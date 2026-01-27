
from libsbml import readSBML
import networkx as nx
import pandas as pd
from xml.dom import minidom
import itertools
import os
import glob

def buildDG(sbml):
    '''
    Usage: reads SBML file, parses reaction and product list
    Returns: networkx directed graph
    '''

    DG = nx.DiGraph()

    document = readSBML(sbml)
    model = document.getModel()

    rxn = (model.getListOfReactions())

    for r in rxn:
        react = [i.getSpecies() for i in r.getListOfReactants()]
        prod = [j.getSpecies() for j in r.getListOfProducts()]
        for r in react:
            for p in prod:
                DG.add_edge(r,p)

    return DG

def getSeedSet(DG, maxComponentSize=5):
    '''
    Usage: takes input networkX directed graph
    Returns: SeedSet dictionary{seedset:confidence score}
    Implementation follows literature description, 
    Improves upon NetCooperate module implementation which erroneously discards certain cases of SCCs (where a smaller potential SCC lies within a larger SCC)
    '''

    SCC = nx.strongly_connected_components(DG)
    SeedSetConfidence = dict()

    for cc in SCC:
        cc_temp = list(cc)

        if len(cc_temp) > int(maxComponentSize):
            continue

        elif len(cc_temp) == 1:
            if DG.in_degree(cc_temp[0]) == 0:
                SeedSetConfidence[cc_temp[0]] = 1.0    

        else:
            for node in cc_temp:
                for edge in DG.in_edges(node):
                    if edge[0] not in cc_temp:
                        cc_temp = []
            for node in cc_temp:
                SeedSetConfidence[node] = 1 / len(cc_temp) 
    
    SeedSet = SeedSetConfidence.keys()
    nonSeedSet = list(set(DG.nodes()) - set(SeedSet))

    return SeedSetConfidence, SeedSet, nonSeedSet

def get_ptm(A_conf, B_nonseedset):
    '''
    Usage: Get potentially transferable metabolites from donor to receptor; A: Receptor, B: Donor
    Returns: [ptm_bigg_ids, ...] (list)
    '''

    SeedA = set(A_conf.keys())
    nonSeedB = set(B_nonseedset)
    intersect_seedA_nonseedB = SeedA.intersection(nonSeedB)

    return list(intersect_seedA_nonseedB)

xml_files = glob.glob('./*.xml') # List of XML files

results = []
for pair in itertools.combinations(xml_files, 2):
    A_filepath, B_filepath = pair
    A_graph, B_graph = buildDG(A_filepath), buildDG(B_filepath)
    A_model_name = minidom.parse(A_filepath).documentElement.getElementsByTagName("model")[0].getAttribute("metaid")
    B_model_name = minidom.parse(B_filepath).documentElement.getElementsByTagName("model")[0].getAttribute("metaid")
    A_conf, A_seedset, A_nonseedset = getSeedSet(A_graph, 5)
    B_conf, B_seedset, B_nonseedset = getSeedSet(B_graph, 5)

    ptm_a_to_b = get_ptm(B_conf, A_nonseedset)
    ptm_b_to_a = get_ptm(A_conf, B_nonseedset)

    bigg_db = pd.read_csv('./bigg_models_metabolites.txt', sep='\t')
    bigg_dict = pd.Series(bigg_db.set_index('bigg_id').T.to_dict('list'))

    for name in ptm_a_to_b:
        name = name[2:]
        results.append([A_model_name, B_model_name, name, bigg_dict[name][0], bigg_dict[name][1], bigg_dict[name][2], bigg_dict[name][3]])

    for name in ptm_b_to_a:
        name = name[2:]
        results.append([B_model_name, A_model_name, name, bigg_dict[name][0], bigg_dict[name][1], bigg_dict[name][2], bigg_dict[name][3]])

with open("./wwtp955_shared_PhyloMInt_PTM.txt", "w") as out:
    out.write('Donor\tReceptor\tbigg_id\tuniversal_bigg_id\tmetabolite_name\tbigg_model_list\tdatabase_links\n')
    for row in results:
        out.write('\t'.join(map(str, row)) + '\n')

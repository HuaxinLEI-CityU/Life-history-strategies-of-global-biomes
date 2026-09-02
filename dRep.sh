#!/bin/bash

dRep dereplicate drep_merine -g all_marine_bins/*fa --genomeInfo marine_checkm_out.csv -p 32
dRep dereplicate drep_surface -g all_surface_bins/*fa --genomeInfo surface_checkm_out.csv -p 32
dRep dereplicate drep_soil -g all_soil_bins/*fa --genomeInfo soil_checkm_out.csv -p 32
dRep dereplicate drep_AS -g all_AS_bins/*fa --genomeInfo AS_checkm_out.csv -p 32

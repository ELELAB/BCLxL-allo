#!/bin/bash

in_pdb=bclxl_puma_3.5nm.pdb
out_pdb=bclxl_puma_3.5nm_processed.pdb

python3.7 renumber_chain.py -i $in_pdb -o $out_pdb --chain B --new-start 130

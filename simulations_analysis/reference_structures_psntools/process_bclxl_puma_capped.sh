#!/bin/bash

in_pdb=bclxl_puma_capped.pdb
out_pdb=bclxl_puma_capped_processed.pdb

python3.7 renumber_chain.py -i $in_pdb -o $out_pdb --chain B --new-start 130

#!/bin/bash

in_pdb=free_2nm_tips3p.pdb
out_pdb=free_2nm_tips3p_processed.pdb

python3.7 assign_chain.py -i $in_pdb -o $out_pdb --chain A --start 1 --end 2642

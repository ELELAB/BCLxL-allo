#!/bin/bash

in_pdb=free_cap_tip3p.pdb
out_pdb=free_cap_tip3p_processed.pdb

python3.7 assign_chain.py -i $in_pdb -o $out_pdb --chain A --start 1 --end 2643

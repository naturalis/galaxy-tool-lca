#!/bin/bash
#if [ $6 == "best_hit" ]
#then
#    lca.py -i $1 -o $2 -b $3 -id $4 -cov $5 -t $6 -tid $7 -tcov $8 -fh $9 -flh "${10}" -minbit "${11}" -fs "${12}"
#fi
#if [ $6 == "only_lca" ]
#then
#    lca.py -i $1 -o $2 -b $3 -id $4 -cov $5 -t $6 -fh $7 -flh $8 -minbit $9 -fs "${10}"
#fi

lca.py -i $1 -o $2 -b $3 -id $4 -cov $5 -t $6 -tid $7 -tcov $8 -fh $9 -flh "${10}" -minbit "${11}" -fs "${12}"

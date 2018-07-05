#!/bin/bash
if [ $6 == "yes" ]
then
    lca.py -i $1 -o $2 -b $3 -id $4 -cov $5 -t $6 -tid $7 -tcov $8 -fh $9 -flh "${10}" -minbit "${11}"
fi
if [ $6 == "no" ]
then
    lca.py -i $1 -o $2 -b $3 -id $4 -cov $5 -t $6 -fh $7 -flh $8 -minbit $9
fi


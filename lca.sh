#!/bin/bash
SCRIPTDIR=$(dirname "$(readlink -f "$0")")
python $SCRIPTDIR"/lca.py" -i $1 -o $2 -b $3 -id $4 -cov $5 -t $6 -tid $7 -tcov $8 -fh $9 -flh "${10}" -minbit "${11}" -fs "${12}"

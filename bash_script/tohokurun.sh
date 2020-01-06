#!/bin/bash
# this script compiles and run all types of the detector designs.
# The result root file will be renamed and saved to result/posrun/TYPE(num)/ directory.

SOURCEDIR=`dirname "$BASH_SOURCE"`/..
echo "
/////////////////////////////////////////////////
// Bash Script tohokurun.sh Starts!
/////////////////////////////////////////////////
"

cd $SOURCEDIR/build
make clean
TYPENAME="TYPE8"
cmake $SOURCEDIR -DDESIGN=$TYPENAME
make -j16
rm -rf $SOURCEDIR/result/tohokurun/
mkdir $SOURCEDIR/result/tohokurun/
cd $SOURCEDIR/result/tohokurun/
$SOURCEDIR/build/simulateAC -m $SOURCEDIR/build/macro/tohokurun/tohokurun_main -t 16 > /dev/null

cd $SOURCEDIR

echo "
////////////////////////////////////////////////
// Bash Script tohokurun.sh finished!
////////////////////////////////////////////////
"

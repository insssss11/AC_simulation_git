#!/bin/bash
# this script compiles and run simulation for aerogel mean free path measurement.
# The result root file will be renamed and saved to result/testbench/ directory.

if [ $# -ne 2 ];then
  echo "Usage : ./testbenchrun.sh [nthread] [condition(1, 2, 3, 4)]"
fi

nthread="$1"
condnum="$2"

SOURCEDIR=`dirname "$BASH_SOURCE"`/..
BUILDDIR="build_testbench$condnum"
echo "
/////////////////////////////////////////////////
// Bash Script testbenchrun.sh Starts!
/////////////////////////////////////////////////
"

# setup
# 1 : led-aerogel-aerogel
# 2 : led-empty-aerogel
# 3 : led-aerogel-empty
# 4 : led-empty-empty

rm -rf $SOURCEDIR/result/testbench/
mkdir $SOURCEDIR/result/testbench/
mkdir $SOURCEDIR/$BUILDDIR
cd $SOURCEDIR/$BUILDDIR
cmake $SOURCEDIR "-DDESIGN=TESTBENCH$condnum"
make -j8
cd $SOURCEDIR/result/testbench/
$SOURCEDIR/$BUILDDIR/simulateAC -m $SOURCEDIR/$BUILDDIR/macro/testbench/testbench_main -t $nthread > /dev/null
mv testbench_355nm.root "setup$condnum.root"

cd $SOURCEDIR
rm -rf $SOURCEDIR/$BUILDDIR

echo "
////////////////////////////////////////////////
// Bash Script testbenchrun.sh finished!
////////////////////////////////////////////////
"

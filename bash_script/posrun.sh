#!/bin/bash
# this script compiles and run given type of the detector designs.
# The result root file will be renamed and saved to result/posrun/TYPE(num)/ directory.

if [ $# -ne 2 ];then
  echo "usage : ./posrun.sh [nthread] [type number]" 
  exit 0
fi

nthread="$1"
TYPENAME="TYPE$2"
cd `dirname "$BASH_SOURCE"`/..
SOURCEDIR=`pwd`
BUILDDIR="build$2"
echo $SOURCEDIR
echo $BUILDDIR
mkdir $SOURCEDIR/$BUILDDIR
cd $SOURCEDIR/$BUILDDIR
g++ -o ProcessTree ../root_macro/ProcessTree.cc `root-config --libs --cflags`

echo "
/////////////////////////////////////////////////
// Bash Script posrun.sh Starts!
/////////////////////////////////////////////////
"

cd $SOURCEDIR/$BUILDDIR
cmake $SOURCEDIR -DDESIGN=$TYPENAME
make
rm -rf $SOURCEDIR/result/posrun/$TYPENAME/
mkdir $SOURCEDIR/result/posrun/$TYPENAME/
cd $SOURCEDIR/result/posrun/$TYPENAME/
$SOURCEDIR/$BUILDDIR/simulateAC -m $SOURCEDIR/$BUILDDIR/macro/posrun/posrun_main -t $nthread > /dev/null
# sort root data file
for rootfile in $(ls *.root)
do
  echo "-----------------sorting $rootfile --------------------"
  $SOURCEDIR/$BUILDDIR/ProcessTree $rootfile
done
cd $SOURCEDIR
echo "---------------deleting build directory....----------"
rm -r $SOURCEDIR/$BUILDDIR

echo "
////////////////////////////////////////////////
// Bash Script posrun.sh finished!
////////////////////////////////////////////////
"

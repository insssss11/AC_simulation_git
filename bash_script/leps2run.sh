#!/bin/bash
# this script compiles and run simulation for Geant4.
# The result root file will be renamed and saved to result/leps2run/ directory.

if [$# -ne 1 ];then
  echo "Usage : ./leps2run.sh [nthread]"
  exit 0;
fi

nthread='$1'
SOURCEDIR=`dirname "$BASH_SOURCE"`/..
BUILDDIR="build_leps2run"
echo "
/////////////////////////////////////////////////
// Bash Script leps2run.sh Starts!
/////////////////////////////////////////////////
"

cd $SOURCEDIR
mkdir $SOURCEDIR/$BUILDDIR
cd $SOURCEDIR/$BUILDDIR
cmake $SOURCEDIR -DDESIGN=LEPS2
make -j8
rm -rf $SOURCEDIR/result/leps2run/
mkdir $SOURCEDIR/result/leps2run/
cd $SOURCEDIR/result/leps2run/
$SOURCEDIR/$BUILDDIR/simulateAC -m $SOURCEDIR/$BUILDDIR/macro/leps2run/leps2run_main -t $nthread > /dev/null

cd $SOURCEDIR
rm -rf $SOURCEDIR/$BUILDDIR

echo "
////////////////////////////////////////////////
// Bash Script leps2run.sh finished!
////////////////////////////////////////////////
"

#! /bin/bash
echo "run this script from bin folder if cmake install was used"
SCRIPT_PATH=`realpath "$0"`
SCRIPT_DIR=`dirname "$SCRIPT_PATH"`
echo in script
echo $SCRIPT_DIR;

export PATH=$SCRIPT_DIR/lib:$PATH
export LD_LIBRARY_PATH=$SCRIPT_DIR/lib:$SCRIPT_DIR/bin:$LD_LIBRARY_PATH
echo appended path, this path is:
echo $PATH
echo appended ld_library_path, this ld... is:
echo $LD_LIBRARY_PATH
$SCRIPT_DIR/bin/decomposer $@

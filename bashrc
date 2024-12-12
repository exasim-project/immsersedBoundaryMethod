
SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

source /opt/OpenFOAM/OpenFOAM-v2306/etc/bashrc

export FOAM_USER_APPBIN=$SOURCE_DIR/code/bin
export FOAM_USER_LIBBIN=$SOURCE_DIR/code/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FOAM_USER_LIBBIN

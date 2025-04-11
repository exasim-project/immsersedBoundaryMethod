
SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#source /opt/OpenFOAM/OpenFOAM-v2306/etc/bashrc

if [ -z "$GIT_ROOT_DIR" ]; then 
    export GIT_ROOT_DIR=$SOURCE_DIR

    #remove $FOAM_USER_LIBBIN from LD_LIBRARY_PATH
    LD_LIBRARY_PATH=$(echo $LD_LIBRARY_PATH | tr ':' '\n' | grep -v -i "$FOAM_USER_LIBBIN" | tr '\n' ':')
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH%:}
    
    #Overwrite OpenFOAM variables
    export FOAM_USER_LIBBIN=$SOURCE_DIR/code/lib
    export FOAM_USER_APPBIN=$SOURCE_DIR/code/bin

    #append the new FOAM_USER_LIBBIN to LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH
else
    echo "GIT_ROOT_DIR already set:  "$GIT_ROOT_DIR
fi

SOURCE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


# Check if OpenFOAMv-2106 has been sourced
if [ "$WM_PROJECT_VERSION" != "v2106" ]; then
    echo "Error: OpenFOAM-v2106 is not sourced. Please source OpenFOAM-v2106 and try again."
    return;
fi


#set GIT_ROOT_DIR if it was not already set
if [ -z "$GIT_ROOT_DIR" ]; then 
    export GIT_ROOT_DIR=$SOURCE_DIR
    
    #Overwrite OpenFOAM variables
    export FOAM_USER_LIBBIN=$SOURCE_DIR/code/lib
    export FOAM_USER_APPBIN=$SOURCE_DIR/code/bin
    export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$LD_LIBRARY_PATH

    export MMF_FOAM_DIR=$SOURCE_DIR/code/multiMultiFoam
    export TPF_FOAM_DIR=$SOURCE_DIR/code/twoPhaseFlow/code/ThirdParty
else
    echo "GIT_ROOT_DIR already set:  "$GIT_ROOT_DIR
    echo "MMF_FOAM_DIR already set:  "$MMF_FOAM_DIR
    echo "TPF_FOAM_DIR already set:  "$TPF_FOAM_DIR
fi


#set LOCAL_SCRIPT_PATH if it was not already set
if [ -z "$LOCAL_SCRIPT_PATH" ]; then 
    if [ ! -e "$SOURCE_DIR/code/mmf_scripts/README.md" ]; then 
        git submodule update --init code/mmf_scripts
    fi
    export LOCAL_SCRIPT_PATH=$SOURCE_DIR/code/mmf_scripts
else
    echo "LOCAL_SCRIPT_PATH already set:  "$LOCAL_SCRIPT_PATH
fi

export PATH=$FOAM_USER_APPBIN:$PATH

echo
echo "Use 'unsetVariables' to clear GIT_ROOT_DIR, MMF_FOAM_DIR, TPF_FOAM_DIR, LOCAL_PYTHON_PATH and LOCAL_SCRIPT_PATH."

unsetVariables(){
    echo "Clearing variables set by $SOURCE_DIR/bashrc"
    unset GIT_ROOT_DIR
    unset MMF_FOAM_DIR
    unset LOCAL_SCRIPT_PATH
    unset TPF_FOAM_DIR
}

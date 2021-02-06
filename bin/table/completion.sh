#/usr/bin/env bash

# TODO: GNU Copyright ...
# Original Author:
# Pedram Ashofteh Ardakani <pedramardakani@pm.me>
# Contributing authors:
# Mohammad Akhlaghi <mohammad@akhlaghi.org>

# A thought on short options. I thing they should not be covered by the
# autocompletion. Because only the advanced users may use them. And it is
# possible to mix them up. So, only those will use the short options who
# know what they are doing. Hence, they will not need the autocompletion
# feature binded to the short options.  However, the short options are
# taken into consideration for suggesting the upcoming commands.

# A though on reporing errors, e.g. invalid filenames. The autocompletion
# feature should not suggest anything in case there is an error in the
# commandline. No errors or messages should be shown as the program in
# charge will handle that. This autocompletion feature is only here to help
# with preventing unintentional mistakes. So, in case there is an invalid
# command on the current commandline, there should be no completion
# suggestions.

# TIP: Run the command below to initialize the bash completion feature for
# this specific program (i.e. astcosmiccal):
# $ source astcosmiccal-completion.bash

# GLOBAL VARIABLES

PREFIX="/usr/local/bin";
ASTFITS="$PREFIX/astfits";
ASTTABLE="$PREFIX/asttable";
db=0 # Set 0 for printing debug messages, else set to 1

# Use extended globs in the case statements if needed
# https://mywiki.wooledge.org/BashGuide/Patterns#Extended_Globs
# shopt -s extglob

#  astquery gaia --dataset=edr3 --center=24,25 --radius=0.1 --output=gaia.fits --column=ra,dec,parallax --quiet -i | awk '/[0-9]+/ {print $2}'

_gnuastro_autocomplete_get_fits_hdu(){
    # Accepts a fits filename as input and echoes its headers
    if [ -f "$1" ]; then
        echo "$($ASTFITS --quiet $1 | awk '{print $2}')"
    fi
}

_gnuastro_autocomplete_list_fits_hdu(){
    # Checks for the current fits file and puts its headers into
    # completion suggestions
    if [ -f "$1"  ]; then
        local list="$(_gnuastro_autocomplete_get_fits_hdu $1)"
        COMPREPLY=($(compgen -W "${list[@]}"))
    fi
}

_gnuastro_autocomplete_list_fits_names(){
    # 'Append' all 'FITS' files in current directory to suggestions. Case
    # insensitive.  The -X option and its filter pattern are explained on
    # bash programmable completion info page: $ info bash programmable
    # 'Appending' seems a good idea because the program might accept
    # multiple input types.

    # For example the 'asttable' program can either accept a fits file or
    # various short/long options as its first argument. In this case,
    # autocompletion suggests both.

    # The completion can not suggest filenames that contain white space in
    # them for the time being.
    COMPREPLY+=($(compgen -f -X "!*.[fF][iI][tT][sS]" -- "$word"))
}

_gnuastro_autocomplete_expect_number(){
    # Prompt the user that this option only accepts a number
    echo "Pass"
}

_gnuastro_autocomplete_get_fits_name(){
    # Get the first fits file among the command line and put it into the
    # $comp_fits_name variable
    # TODO: How about all other fits file extensions?
    local file_name="$(echo ${COMP_WORDS[@]} | awk -v regex="[a-zA-Z0-9]*.[fF][iI][tT][sS]" 'match($0, regex) {print substr($0, RSTART, RLENGTH)}')"
    if [ -f "$file_name" ]; then
        # Check if file_name is actually an existing fits file. This
        # prevents other functions from failing and producing obscure error
        # messages
        echo "$file_name"
        # Note that should not be an 'else' statement with 'exit' error
        # code. Because this function is checking the presence of a fits
        # file everytime bash completion is provoked. Then it will return
        # error if there is no fits name and break functionality.
    fi
}

_gnuastro_autocomplete_get_fits_columns(){
    # Checks if the argument contains a valid file. Does not check for its
    # extension. Then, reads the column names using the asttable program
    # and echoes the resulting STR.
    if [ -f "$1" ]; then
        # Force 'awk' to read after the second line of 'asttable' output,
        # because the second line contains the filename. The filename might
        # start with numbers. If so, there will be an unwanted '(hdu:'
        # printed in the results. Here, 'awk' will print the second column
        # in lines that start with a number.
        echo "$($ASTTABLE --information $1 | awk 'NR>2' | awk '/^[0-9]/ {print $2}')"
    fi
}

_gnuastro_autocomplete_list_fits_columns(){
    # Accept a fits file name as the first argument ($1). Read and suggest
    # its column names. If the file does not exist, pass.
    if [ -f "$1" ]; then
        local list="$(_gnuastro_autocomplete_get_fits_columns $1)"
        COMPREPLY=($(compgen -W "$list"))
    fi
}

_gnuastro_autocomplete_get_file(){
    # The last file name (.txt/.fits) in COMP_LINE
    echo "Pass"
}

# just find the short commands
# astconvolve --help | awk -v pattern="^ *-([a-z]|[A-Z])" 'match($0, pattern) {print $0}'

_gnuastro_autocomplete_list_options(){
    # Accept the command name and its absolute path, run the --help option
    # and 'append' all long options to the current suggestions. 'Appending'
    # seems a good idea because the program might accept multiple input
    # types. For example the 'asttable' program can either accept a fits
    # file or various short/long options as its first argument. In this
    # case, autocompletion suggests both.
    local list=$("$1" --help | awk -v regex=" --+[a-zA-Z0-9]*=?" 'match($0, regex) {print substr($0, RSTART, RLENGTH)}')
    COMPREPLY+=($(compgen -W "$list" -- "$word"))
}

_gnuastro_asttable_completions(){

    # TODO: @@
    local PROG_NAME="asttable";

    local PROG_ADDRESS="$PREFIX/$PROG_NAME";

    # Initialize the completion response with null
    COMPREPLY=();

    # Variable "word", is the current word being completed
    local word="${COMP_WORDS[COMP_CWORD]}";

    # Variable "prev" is the word just before the current word
    local prev="${COMP_WORDS[COMP_CWORD-1]}";

    # A quick check to see if there is already a fits file name invoked in
    # the current commandline. This means the order of commands does matter
    # in this bash completion. If we do not want this, we should implement
    # another method for suggesting completions.
    local fits_name="$(_gnuastro_autocomplete_get_fits_name)"

    # TODO: Prettify the code syntax, shorter ones on top
    case "$prev" in
        asttable)
            _gnuastro_autocomplete_list_fits_names
            _gnuastro_autocomplete_list_options $PROG_ADDRESS
            ;;
        -i|--information|-w|--wcsfile)
            if [ -f "$fits_name" ]; then
                # The user has entered a valid fits file name. So keep on
                # with suggesting all other options at hand.
                _gnuastro_autocomplete_list_options $PROG_ADDRESS
            else
                # Check if the user has already specified a fits file. If
                # the _gnuastro_autocomplete_get_file_name echoes an empty
                # response, it means no fits files were specified.
                _gnuastro_autocomplete_list_fits_names
            fi
            ;;
        -c|--column|-r|--range|-s|--sort)
            # The function below, checks if the user has specified a fits
            # file in the current commandline. If not, there will be no
            # response from autocompletion. This might alert the user that
            # something is going wrong.
            _gnuastro_autocomplete_list_fits_columns "$fits_name"
            ;;
        -W|--wcshdu)
            # Description is same as the '--column' option.
            _gnuastro_autocomplete_list_fits_hdu "$fits_name"
            ;;
        -b|--noblank) ;;
        -h|--hdu) ;;
        *) _gnuastro_autocomplete_list_options $PROG_ADDRESS ;;
    esac

    if [[ ! "${COMPREPLY[@]}" =~ "=" ]]; then
        # '[[' and 'compopt' work for bash 4+
        # TODO: Find workaround for older bash
        compopt +o nospace
    fi

    # Debugging purpose:
    if [ $db -eq 0 ]; then
        cat <<EOF

*** DEBUG ***
>>> prev: '$prev' -- \$3: '$3'
>>> word: '$word' -- \$2: '$2'
>>> fits_name: '$fits_name'
EOF
        printf ">>> line: ${COMP_LINE[@]}"
    fi

}

complete -F _gnuastro_asttable_completions -o nospace asttable

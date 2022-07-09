#!/bin/bash
shopt -s extglob

path=${1%%*([0-9a-zA-Z-._])}  # remove filename
path=${path//[ ]/' '}         # enclose space in quotes
mode=$2

# Debug
# echo $path
# echo $2

if [ $mode = "pull" ]; then
    git -C "$path" pull
elif [ $mode = "push" ]; then
    git -C "$path" commit -a -m "Committed by TeXstudio"
    git -C "$path" push -u origin "main"
elif [ $mode = "stash" ]; then
    git -C "$path" stash
fi

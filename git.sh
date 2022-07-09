#!/bin/bash
shopt -s extglob

path=${1%%*([0-9a-zA-Z-._])}
path=${path//[ ]/' '}
mode=$2

echo $path
echo $2

if [ $mode = "pull" ]
then
    git -C "$path" pull
else
    git -C "$path" commit -a -m "Committed by TeXstudio"
    git -C "$path" push -u origin "main"
fi

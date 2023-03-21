#!/bin/bash

# This script scans all subdirectories for files containing "sync-conflict" in the name and deletes them.

# Find all files containing "sync-conflict" in the name and delete them
find . -type f -name "*sync-conflict*" -exec rm {} \;

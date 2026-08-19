#!/bin/bash

echo "In space, no one can hear you scream."

for file in *\ *; do
    [ -e "$file" ] || continue
    mv -- "$file" "${file// /_}"
done

#!/bin/bash

# Submit all scripts ending in W.sh
for file in *W.sh; do
    if [[ -f "$file" ]]; then
        echo "Submitting $file"
        sbatch "$file"
    fi
done

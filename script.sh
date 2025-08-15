#!/bin/bash
for i in {1..10000};
do
    python3 make_input.py 100000
    echo -n "*"
    output=$(./tripod_exec adjacencies.txt simplices.txt triangles.txt mybasename)
    echo -n -e "\r"$i
    if echo "$output" | grep -q "error|unsuccessful";
    then
        echo "error found"
        break
    fi
done
echo

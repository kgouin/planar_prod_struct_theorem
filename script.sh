#!/bin/bash
for i in {1..10};
do
    python3 make_input.py 1000
    echo -n "*"
    output=$(./tripod_exec)
    echo -n -e "\r"$i
    if echo "$output" | grep -q "error|unsuccessful";
    then
        echo "error found"
        break
    fi
done
echo

#!/bin/bash
clear
make
./tripod_exec > new_output.txt
diff original_output.txt new_output.txt

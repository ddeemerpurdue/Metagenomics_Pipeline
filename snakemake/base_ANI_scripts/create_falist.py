import os
import sys

line = []
for file in os.listdir(sys.argv[1]):
    if file.endswith('.fa'):
        line.append(file)
line = sorted(line)
with open(sys.argv[2], 'w') as o:
    for l in line:
        o.write(str(sys.argv[1]) + str(l))
        o.write('\n')


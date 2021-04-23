import sys
import numpy as np

file_name = sys.argv[1]
size = int(sys.argv[2])

file = open(file_name, "w")
newline = "\n"
for i in range(size):
    for j in range(size):
        a = np.random.uniform(-5.0, 5.0)
        s = '%.12f'%a
        file.write(s)
        file.write(" ")
    file.write(newline)

file.close()

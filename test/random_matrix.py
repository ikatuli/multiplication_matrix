import random
import numpy

maximum=30

random.seed();

m=random.randint(1, maximum)
n=random.randint(1, maximum)
a=numpy.array([[random.randint(1, maximum)for i in range(n)] for i in range(m)])

file = open("A.txt", "w+")
file.write(str(m)+" "+str(n)+'\n')
for i in a:
    for j in i:
        file.write(str(j))
        file.write(" ")
    file.write("\n")

m=random.randint(1, maximum)
b=numpy.array([[random.randint(1, maximum)for i in range(m)] for i in range(n)])

file = open("B.txt", "w+")
file.write(str(n)+" "+str(m)+'\n')
for i in b:
    for j in i:
        file.write(str(j))
        file.write(" ")
    file.write("\n")

C = a.dot(b)

file = open("C.txt", "w+")
file.write(str(len(C[0]))+" "+str(len(C))+'\n')
for i in C:
    for j in i:
        file.write(str(j))
        file.write(" ")
    file.write("\n")

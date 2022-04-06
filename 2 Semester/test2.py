# import re
import numpy as np
import math

with open("test.srt") as f:
    d = f.readlines()

x = " ".join(d)
print(x)
words = []
oldwords = x.split(" ")
newwords = []

for word in oldwords:
    temp = list(word)
    y = list(np.random.randint(0, len(temp), math.ceil(len(temp) / 4)))
    y.sort(reverse=True)
    for i in y:
        temp[i] = ""
    words = "".join(temp)
    newwords.append(words)

print(newwords)
# y = np.random.randint(0, len(x), round(len(x)/3))
# new = list(x)
z = []
# print(neww)

# for i in y:
#     if new[i] != ' ':
#         z.append(i)
#
# for i in z:
#     new[i] = ""

print(" ".join(newwords))

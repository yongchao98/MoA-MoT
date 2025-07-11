from re import *
print(*(sorted(list(set(map(int, findall(r'\d+', input()))))) or ["NO"]))
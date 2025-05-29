import math

# Side lengths
a = 5
b = 6
c = 7

# Semi-perimeter
s = (a + b + c) / 2

# Area using Heron's formula
A = math.sqrt(s * (s - a) * (s - b) * (s - c))

# Inradius
r = A / s

print(r)
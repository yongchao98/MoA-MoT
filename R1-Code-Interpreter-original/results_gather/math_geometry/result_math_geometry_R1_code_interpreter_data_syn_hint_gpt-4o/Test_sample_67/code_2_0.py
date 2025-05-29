import math

# Given lengths
YX = 60
XZ = 80

# Calculate the hypotenuse YZ using the Pythagorean theorem
YZ = math.sqrt(YX**2 + XZ**2)

# Calculate the coordinates of X
a = (YZ**2 + YX**2 - XZ**2) / (2 * YZ)
b = math.sqrt(YX**2 - a**2)

# Calculate the coordinates of Z
c = YZ

# Calculate WZ
WZ = abs(c - a)

# Output the result
print(WZ)
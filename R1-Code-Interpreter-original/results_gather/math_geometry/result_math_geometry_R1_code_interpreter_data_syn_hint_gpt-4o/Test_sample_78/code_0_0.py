import math

# Given values
PQ = 12
cos_Q = 0.4

# Calculate PR
PR = PQ / cos_Q

# Calculate QR using the Pythagorean theorem
QR = math.sqrt(PR**2 - PQ**2)

print(QR)
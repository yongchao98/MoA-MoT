import math

# Coordinates of C and C'
C = (-3, 2)
C_prime = (3, 2)

# Calculate the distance
distance = math.sqrt((C_prime[0] - C[0])**2 + (C_prime[1] - C[1])**2)
print(distance)
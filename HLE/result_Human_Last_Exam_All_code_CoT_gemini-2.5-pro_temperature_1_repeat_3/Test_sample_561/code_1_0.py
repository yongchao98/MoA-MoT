import math

# N is the number of self-similar copies. In a piano octave, there are 5 black keys,
# so the keyboard is replaced by 5 smaller copies of itself.
N = 5

# r is the scaling factor.
# The original keyboard is 3 units wide.
# There are 7 white keys, so each is 3/7 units wide.
# A black key is half as wide: (1/2) * (3/7) = 3/14 units.
# The scaling factor r = (new width) / (original width) = (3/14) / 3 = 1/14.
# The inverse of the scaling factor, 1/r, is therefore 14.
inv_r = 14

# The Minkowski–Bouligand dimension D for a self-similar fractal is log(N) / log(1/r).
dimension = math.log(N) / math.log(inv_r)

print("The Minkowski–Bouligand dimension (D) is calculated using the formula: D = log(N) / log(1/r)")
print(f"Here, N = {N} (the number of black keys).")
print(f"The scaling factor r = 1/14, so 1/r = {inv_r}.")
print("\nThe final equation is:")
print(f"D = log({N}) / log({inv_r})")
print(f"D = {math.log(N)} / {math.log(inv_r)}")
print(f"D = {dimension}")
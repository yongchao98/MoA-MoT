import math

# Step 1: Define the number of self-similar copies (N)
# An octave has 5 black keys.
N = 5

# Step 2: Define the scaling factor (s)
# Keyboard dimensions: width=3, height=1
# White key width = 3 / 7
# Black key width = (1/2) * (3/7) = 3/14
# Black key height = 9/14 * 1 = 9/14
# A scaled keyboard (3r x 1r) must fit inside a black key (3/14 x 9/14).
# Width: 3r <= 3/14 => r <= 1/14
# Height: 1r <= 9/14 => r <= 9/14
# The limiting scaling ratio is r = 1/14.
# The scaling factor s is the reciprocal of r.
s = 14

# Step 3: Calculate the Minkowski-Bouligand dimension (D)
# D = log(N) / log(s)
# The base of the logarithm does not matter as long as it's consistent.
D = math.log(N) / math.log(s)

# Step 4: Print the explanation and the result
print("The Minkowski-Bouligand dimension (D) is calculated for a self-similar fractal using the formula:")
print("D = log(N) / log(s)")
print("-" * 30)
print(f"Number of self-similar copies (N): {N}")
print(f"Scaling factor (s): {s}")
print(f"The final equation is: D = log({N}) / log({s})")
print(f"The calculated dimension is: {D}")
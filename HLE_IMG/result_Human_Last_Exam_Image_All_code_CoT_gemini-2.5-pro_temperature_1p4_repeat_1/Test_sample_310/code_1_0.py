# Step 1: Define the parameters for the missing simulation.
# The missing plot is for a linear polymer.
# The initial degree of polymerization (N) can be identified from the plots as the location of the initial single peak.
N = 20

# The degree of destruction 'a' for the missing plot was identified as m=4.
# a = m / 25
m = 4
a = m / 25

# The fraction of broken bonds 'p' is equal to the degree of destruction 'a' for random scission.
p = a

# Step 2: Use the formula for the number-average degree of polymerization (Nn).
# The formula is Nn = N / (1 + p * (N - 1))
Nn = N / (1 + p * (N - 1))

# Step 3: Print the calculation and the final result.
print(f"The number-average degree of polymerization (Nn) for the missing plot is calculated as follows:")
print(f"Nn = N / (1 + a * (N - 1))")
print(f"Nn = {N} / (1 + {a} * ({N} - 1))")
print(f"Nn = {N} / (1 + {a} * {N-1})")
print(f"Nn = {N} / (1 + {a * (N-1)})")
print(f"Nn = {N} / {1 + a * (N-1)}")
print(f"Nn = {Nn}")
# Parameters for the missing simulation identified from the analysis.
# N is the initial degree of polymerization.
N = 20
# The missing plot corresponds to a linear polymer simulation with m=4.
m = 4

# Calculate the degree of destruction, a.
a = m / 25.0

# Calculate the number-average degree of polymerization, N_n, using the standard formula
# for random scission of a monodisperse linear polymer.
N_n = N / (1 + a * (N - 1))

# Print the equation with all the numbers substituted for clarity.
print("Calculation for the number-average degree of polymerization (N_n):")
print(f"Formula: N_n = N / (1 + a * (N - 1))")
print(f"Substituting the values for the missing plot:")
print(f"N = {N}")
print(f"a = {m}/25 = {a}")
print(f"N_n = {N} / (1 + {a} * ({N} - 1))")
print(f"N_n = {N} / (1 + {a * (N - 1)})")
print(f"N_n = {N} / {1 + a * (N - 1)}")
print(f"The final result is: {N_n}")
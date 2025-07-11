# Parameters for the polymer simulation
# N is the initial degree of polymerization. From the plots, the initial polymer has a size of 25.
N = 25

# m is the integer that defines the degree of destruction, a = m/N.
# Based on the analysis of the provided plots, the missing simulation corresponds to m=4 for a linear polymer.
m = 4

# Calculate the degree of destruction 'a'
a = m / N

# Calculate the number of scissible bonds in the initial polymer chain
N_minus_1 = N - 1

# Calculate the number-average degree of polymerization (Nn) using the formula for random scission of linear polymers.
# The formula is Nn = N / (1 + a * (N - 1))
Nn = N / (1 + a * N_minus_1)

# Print the calculation step-by-step as requested.
print(f"The number-average degree of polymerization (Nn) for the missing plot is calculated as follows:")
print(f"The formula for a linear polymer is: Nn = N / (1 + a * (N - 1))")
print(f"The initial degree of polymerization is N = {N}.")
print(f"The missing plot corresponds to m = {m}.")
print(f"The degree of destruction is a = m / N = {m} / {N} = {a}.")
print(f"Plugging the values into the formula:")
print(f"Nn = {N} / (1 + {a} * ({N} - 1))")
print(f"Nn = {N} / (1 + {a} * {N_minus_1})")
print(f"Nn = {N} / (1 + {a * N_minus_1})")
print(f"Nn = {N} / {1 + a * N_minus_1}")
print(f"Nn = {Nn}")
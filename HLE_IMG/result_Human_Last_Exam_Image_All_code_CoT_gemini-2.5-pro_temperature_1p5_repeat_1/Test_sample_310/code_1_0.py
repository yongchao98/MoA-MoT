# Parameters for the simulation
# Initial degree of polymerization (deduced from plots and problem description)
N = 25
# Parameter m for the degree of destruction 'a' for the missing plot
# My analysis concluded the missing plot is for the linear polymer with the lowest
# degree of destruction, which corresponds to m=1.
m = 1

# Calculate the degree of destruction 'a'
a = m / 25.0

# Define the formula for the number-average degree of polymerization for a linear polymer
# N_n = N / (1 + (N-1)*a)
numerator = N
denominator = 1 + (N - 1) * a

# Calculate the final value
N_n = numerator / denominator

# Output the equation with the values plugged in, and the final result
print(f"The missing simulation corresponds to a linear polymer with a = {m}/25 = {a}.")
print("The number-average degree of polymerization, Nn, is calculated as:")
print(f"Nn = N / (1 + (N - 1) * a)")
print(f"Nn = {N} / (1 + ({N} - 1) * {a})")
print(f"Nn = {N} / (1 + {N-1} * {a})")
print(f"Nn = {N} / (1 + { (N-1) * a })")
print(f"Nn = {N} / {denominator}")
print(f"Nn = {N_n}")
print(f"\nFinal Answer: {N_n}")
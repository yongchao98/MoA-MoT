# Define the pollination rates
self_pollination_rate = 0.05
cross_pollination_rate = 0.95

# Define the fraction of resistant offspring for each pollination type
# For self-pollination (Wt/Ins x Wt/Ins), resistant offspring are Wt/Ins and Ins/Ins (3/4)
resistant_fraction_self = 0.75
# For cross-pollination (Wt/Ins x Wt/Wt), resistant offspring are Wt/Ins (1/2)
resistant_fraction_cross = 0.50

# Calculate the contribution from each pollination type
contribution_from_self = self_pollination_rate * resistant_fraction_self
contribution_from_cross = cross_pollination_rate * resistant_fraction_cross

# Calculate the total theoretical probability of resistant offspring
total_resistant_probability = contribution_from_self + contribution_from_cross

# Convert the probability to a percentage
total_resistant_percentage = total_resistant_probability * 100

# Print the final equation and the result
print("The calculation for the theoretical percentage of drought-resistant offspring is:")
print(f"({self_pollination_rate} * {resistant_fraction_self}) + ({cross_pollination_rate} * {resistant_fraction_cross}) = {total_resistant_probability}")
print(f"Theoretically, {total_resistant_percentage:.2f}% of the offspring should be drought-resistant.")

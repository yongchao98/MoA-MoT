# The temporal horizon constant
temporal_horizon = 48

# A list to store the found temporal fixed points
fixed_points = []

# We are looking for positive integer temporal fixed points k less than 100.
# The condition for a fixed point is |k+k| mod h(T) = |k-k| mod h(T),
# which simplifies to (2 * k) % temporal_horizon == 0.
for k in range(1, 100):
    if (2 * k) % temporal_horizon == 0:
        fixed_points.append(k)

# Calculate the sum of all found fixed points
total_sum = sum(fixed_points)

# Create the equation string by joining the numbers with ' + '
equation_str = " + ".join(map(str, fixed_points))

# Print the final equation with the result
print(f"{equation_str} = {total_sum}")
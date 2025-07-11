# The maximal distance can be found by choosing a specific arrangement of observers
# and a specific movement pattern for the snail, as explained in the plan.

# We model the snail's movement over 21 intervals, each 1/3 of a minute long.
# Let x_i be the distance traveled in the i-th interval.
distances_per_interval = []

# According to the optimal strategy, the snail moves 1 meter only during
# the "even" intervals (2nd, 4th, 6th, etc.) and stays put during the "odd" ones.
for i in range(1, 22):  # i from 1 to 21
    if i % 2 == 0:
        # even interval
        distances_per_interval.append(1)
    else:
        # odd interval
        distances_per_interval.append(0)

# The total maximal distance is the sum of distances from all 21 intervals.
total_distance = sum(distances_per_interval)

# To fulfill the user's request, we will print the full equation.
equation_parts = [str(d) for d in distances_per_interval]
equation_str = " + ".join(equation_parts)

print("The maximal possible distance is found by summing the distances traveled in each 1/3-minute interval:")
print(f"{equation_str} = {total_distance}")

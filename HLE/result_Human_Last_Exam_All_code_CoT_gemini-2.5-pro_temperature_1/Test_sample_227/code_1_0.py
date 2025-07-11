# Each term in the sum represents the distance the snail travels during a specific time interval
# where it is being watched by exactly one observer. According to our strategy, the snail
# covers exactly 1 meter in each of these intervals.

# Distance traveled in [0, 0.5) minutes
d1 = 1
# Distance traveled in (1, 1.5) minutes
d2 = 1
# Distance traveled in [1.5, 2.5) minutes
d3 = 1
# Distance traveled in [2.5, 3.5) minutes
d4 = 1
# Distance traveled in [3.5, 4.5) minutes
d5 = 1
# Distance traveled in [4.5, 5.5) minutes
d6 = 1
# Distance traveled in (5.5, 6) minutes
d7 = 1
# Distance traveled in (6.5, 7] minutes
d8 = 1

# The snail rests during intervals [0.5, 1] and [6, 6.5] where multiple observers are watching.

# The total distance is the sum of the distances from each movement phase.
total_distance = d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8

# Print the equation with each number, as requested.
print(f"{d1} + {d2} + {d3} + {d4} + {d5} + {d6} + {d7} + {d8} = {total_distance}")

# This script determines the four-digit sequence by analyzing the simulations.
# The analysis and resulting codes are determined as follows.

# Parameter encoding map:
# w1 -> 1, w2 -> 2, w3 -> 3
# u1 -> 4, u2 -> 5, u3 -> 6
# h1 -> 7, h2 -> 8, h3 -> 9

# Analysis of Simulation 1:
# - h (intensity) is increasing: h1 < h2 < h3 (Rule 2)
# - u (tilt) is constant: u1 = u2 = u3 (Rule 1)
# - w (spread) has a unique middle value: w1 = w3, w2 is different (Rule 3)
# Unique parameter code for Simulation 1 is for w2.
sim1_code = 2

# Analysis of Simulation 2:
# - h (intensity) is constant: h1 = h2 = h3 (Rule 1)
# - u (tilt) is decreasing: u1 > u2 > u3 (Rule 2)
# - w (spread) has a unique left value: w2 = w3, w1 is different (Rule 3)
# Unique parameter code for Simulation 2 is for w1.
sim2_code = 1

# Analysis of Simulation 3:
# - w (spread) is constant: w1 = w2 = w3 (Rule 1)
# - u (tilt) is increasing: u1 < u2 < u3 (Rule 2)
# - h (intensity) has a unique middle value: h1 = h3, h2 is different (Rule 3)
# Unique parameter code for Simulation 3 is for h2.
sim3_code = 8

# Analysis of Simulation 4:
# - h (intensity) is constant: h1 = h2 = h3 (Rule 1)
# - w (spread) is increasing: w1 < w2 < w3 (Rule 2)
# - u (asymmetry) has a unique right value: u1 = u2, u3 is different (Rule 3)
# Unique parameter code for Simulation 4 is for u3.
sim4_code = 6

# Printing the final sequence
print("The four-digit sequence based on the unique parameter for each simulation is:")
print(f"{sim1_code}{sim2_code}{sim3_code}{sim4_code}")
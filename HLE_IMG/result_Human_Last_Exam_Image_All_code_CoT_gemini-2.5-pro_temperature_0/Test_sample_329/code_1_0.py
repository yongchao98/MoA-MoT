# This script calculates and prints the four-digit sequence based on the analysis of the triple-slit simulations.

# The analysis for each simulation identifies the unique parameter according to Rule 3.
# Rule 3: One parameter differs for only one slit.

# The encoding for each unique parameter is:
# w1: 1, w2: 2, w3: 3
# u1: 4, u2: 5, u3: 6
# h1: 7, h2: 8, h3: 9

# Simulation 1: Unique parameter is w2.
sim1_code = 2

# Simulation 2: Unique parameter is h3.
sim2_code = 9

# Simulation 3: Unique parameter is u2.
sim3_code = 5

# Simulation 4: Unique parameter is w1.
sim4_code = 1

# The final four-digit sequence is formed by the individual codes.
print(f"{sim1_code}{sim2_code}{sim3_code}{sim4_code}")
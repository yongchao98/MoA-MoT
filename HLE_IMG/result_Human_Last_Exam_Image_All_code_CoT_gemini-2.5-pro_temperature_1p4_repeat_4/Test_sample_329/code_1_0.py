# Based on the visual analysis of the four triple-slit simulations, we can deduce the unique parameter for each.
# The analysis for each simulation is as follows:

# Simulation 1:
# - Rule 1 (Constant): Beam direction (u) is constant for all slits.
# - Rule 2 (Monotonic): Brightness (h) decreases from left to right.
# - Rule 3 (Unique): The width (w) of the central slit's diffraction pattern is wider, meaning its slit width (w2) is uniquely narrow.
sim1_code = 2  # Code for w2

# Simulation 2:
# - Rule 1 (Constant): Brightness (h) is constant for all slits.
# - Rule 2 (Monotonic): Beam direction (u) changes monotonically.
# - Rule 3 (Unique): The width (w) of the left slit's diffraction pattern is wider, meaning its slit width (w1) is uniquely narrow.
sim2_code = 1  # Code for w1

# Simulation 3:
# - Rule 1 (Constant): The width (w) of the diffraction patterns is constant.
# - Rule 2 (Monotonic): Beam direction (u) changes monotonically.
# - Rule 3 (Unique): The brightness (h) of the central slit is uniquely high.
sim3_code = 8  # Code for h2

# Simulation 4:
# - Rule 1 (Constant): The width (w) of the diffraction patterns is constant.
# - Rule 2 (Monotonic): Brightness (h) increases from left to right.
# - Rule 3 (Unique): The beam direction (u) of the rightmost slit is unique, as the other two are parallel.
sim4_code = 6  # Code for u3

# The final four-digit sequence is formed by combining the codes.
print(f"{sim1_code}{sim2_code}{sim3_code}{sim4_code}")
# The plan is to determine the unique parameter code for each simulation
# based on the visual evidence and the rules provided.

# Simulation 1:
# - Rule 1 (Identical): Wave number 'u' (parallel beams).
# - Rule 2 (Monotonic): Height 'h' (brightness decreases from left to right).
# - Rule 3 (Unique): Width 'w'. The central beam's spread is widest, meaning its slit 'w2' is narrowest.
sim1_code = 2
print(f"The unique parameter for Simulation 1 is w2, which has the code: {sim1_code}")

# Simulation 2:
# - Rule 1 (Identical): Wave number 'u' (parallel beams).
# - Rule 2 (Monotonic): Width 'w' (beam spread increases left to right, so slit width decreases).
# - Rule 3 (Unique): Height 'h'. The central slit is dimmer than the outer two. 'h2' is unique.
sim2_code = 8
print(f"The unique parameter for Simulation 2 is h2, which has the code: {sim2_code}")

# Simulation 3:
# - Rule 1 (Identical): Height 'h' (slits are equally bright).
# - Rule 2 (Monotonic): Wave number 'u' (beams are converging, so angle changes monotonically).
# - Rule 3 (Unique): Width 'w'. Due to symmetry, the central slit 'w2' is the unique one.
sim3_code = 2
print(f"The unique parameter for Simulation 3 is w2, which has the code: {sim3_code}")

# Simulation 4:
# - Rule 3 (Unique): Wave number 'u'. The left beam 'u1' is angled, while the other two are straight.
sim4_code = 4
print(f"The unique parameter for Simulation 4 is u1, which has the code: {sim4_code}")

# The final four-digit sequence is formed by concatenating the codes.
print("\nThe final four-digit sequence is the combination of these codes.")
print(f"Final sequence = {sim1_code}{sim2_code}{sim3_code}{sim4_code}")
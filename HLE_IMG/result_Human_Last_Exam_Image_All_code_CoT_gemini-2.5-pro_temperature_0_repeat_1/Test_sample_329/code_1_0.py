# Step 1: Define the encoding map for reference.
# w1->1, w2->2, w3->3
# u1->4, u2->5, u3->6
# h1->7, h2->8, h3->9

# Step 2: Analyze each simulation to find the unique parameter and its corresponding code.

# Simulation 1:
# - h is monotonic (intensities increase left to right).
# - w is identical (beam widths look the same).
# - Therefore, u must be unique. The pattern asymmetry suggests u2 is the unique parameter.
sim1_code = 5 # for u2

# Simulation 2:
# - h is identical (intensities look the same).
# - w is monotonic (beam widths increase, so slit widths decrease).
# - Therefore, u must be unique. The pattern shift to the right suggests u1 is unique.
sim2_code = 4 # for u1

# Simulation 3:
# - h is unique (h1=h3, h2 is dimmer). The unique parameter is h2.
# - w appears identical (beam widths look the same).
# - Therefore, u must be monotonic.
sim3_code = 8 # for h2

# Simulation 4:
# - h is identical (intensities look the same).
# - u is monotonic (pattern is strongly tilted).
# - Therefore, w must be unique. The right beam is wider, so w3 is unique.
sim4_code = 3 # for w3

# Step 3: Combine the codes into the final four-digit sequence.
# The problem asks to output each number in the final sequence.
print(f"The unique parameter code for simulation 1 is: {sim1_code}")
print(f"The unique parameter code for simulation 2 is: {sim2_code}")
print(f"The unique parameter code for simulation 3 is: {sim3_code}")
print(f"The unique parameter code for simulation 4 is: {sim4_code}")
print(f"The final four-digit sequence is: {sim1_code}{sim2_code}{sim3_code}{sim4_code}")

# Step 1: Analyze Simulation 1 to find its unique parameter code.
# - Intensity (h) increases from left to right (h1 < h2 < h3). This is a strictly increasing parameter (Rule 2).
# - Beam direction (u) is constant for all three slits (beams are parallel). This is the identical parameter (Rule 1).
# - Therefore, slit width (w) must be the parameter that differs for only one slit (Rule 3).
# - Visually, the left beam's envelope is wider than the other two, which are similar. This means w1 is unique.
# - The code for w1 is 1.
sim1_code = 1

# Step 2: Analyze Simulation 2 to find its unique parameter code.
# - Intensity (h) is the same for all three slits. This is the identical parameter (Rule 1).
# - Beam direction (u) changes monotonically across the slits. This is a strictly changing parameter (Rule 2).
# - Therefore, slit width (w) must be the unique parameter (Rule 3).
# - Visually, the middle beam's envelope is wider than the other two, which are similar. This means w2 is unique.
# - The code for w2 is 2.
sim2_code = 2

# Step 3: Analyze Simulation 3 to find its unique parameter code.
# - Slit width (w) appears constant for all three slits (envelopes are the same width). This is the identical parameter (Rule 1).
# - Beam direction (u) changes monotonically across the slits. This is a strictly changing parameter (Rule 2).
# - Therefore, intensity (h) must be the unique parameter (Rule 3).
# - Visually, the middle slit is brighter than the other two, which are similar. This means h2 is unique.
# - The code for h2 is 8.
sim3_code = 8

# Step 4: Analyze Simulation 4 to find its unique parameter code.
# - Intensity (h) is the same for all three slits. This is the identical parameter (Rule 1).
# - Slit width (w) changes monotonically (the beam envelope gets narrower from left to right). This is a strictly changing parameter (Rule 2).
# - Therefore, beam direction (u) must be the unique parameter (Rule 3).
# - Visually, the first two beams are parallel, while the third is tilted. This means u3 is unique.
# - The code for u3 is 6.
sim4_code = 6

# Step 5: Combine the codes to form the final four-digit sequence.
final_sequence = f"{sim1_code}{sim2_code}{sim3_code}{sim4_code}"

# Step 6: Print the final answer.
print(final_sequence)
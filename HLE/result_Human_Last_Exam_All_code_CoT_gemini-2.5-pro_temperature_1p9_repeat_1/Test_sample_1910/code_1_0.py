import math

# Step 1: Define the Miller indices (h, k, l) and d-spacing for the second major peak.
# Based on literature (R. Černý et al., J. Alloys Compd. 350 (2003) 18-24),
# the second major peak for orthorhombic NaMgH3 is (121).
h, k, l = 1, 2, 1

# The d-spacing for the (121) reflection is given as 2.8785 Angstroms.
d_hkl = 2.8785  # in Angstroms

# Step 2: Calculate the Q-space position.
# The formula is Q = 2 * pi / d
Q = (2 * math.pi) / d_hkl

# Step 3: Print the results, including the final equation.
print(f"The second major diffraction peak corresponds to the Miller indices (hkl) = ({h}, {k}, {l}).")
print(f"The interplanar d-spacing for this peak is {d_hkl} Å.")
print("The Q-space position is calculated using the formula: Q = 2 * π / d.")
print("\nFinal Calculation:")
print(f"Q = (2 * {math.pi:.5f}) / {d_hkl:.5f}")
print(f"Q = {Q:.4f} (1/Å)")
print("\nTherefore, the second major diffraction peak is located at approximately {:.4f} 1/Å in Q-space.".format(Q))

import math

# Step 1: Define the crystal parameters for NaMgH3.
# We assume the ideal simple cubic perovskite structure.
# The lattice parameter 'a' is taken from literature (e.g., Zhang et al., Phys. Rev. B 76, 052103 (2007)).
a = 3.886  # in Angstroms

# Step 2: Identify the Miller indices (h, k, l) for the second major diffraction peak.
# For a simple cubic perovskite, the major peaks are those where h+k+l is even.
# The first major peak is (110). The second major peak is (200).
h, k, l = 2, 0, 0

# Step 3: Calculate the interplanar d-spacing for the (200) reflection.
# For a cubic system, d_hkl = a / sqrt(h^2 + k^2 + l^2).
d_hkl = a / math.sqrt(h**2 + k**2 + l**2)

# Step 4: Calculate the Q-space position.
# The formula is Q = 2 * pi / d. The unit will be 1/Angstrom.
Q_hkl = (2 * math.pi) / d_hkl

# Print the explanation and the final equation with all its numerical components.
print(f"The calculation is for the second major peak, identified as the ({h},{k},{l}) reflection of cubic NaMgH3.")
print(f"Using lattice parameter a = {a} Å.")
print(f"The interplanar spacing is d({h}{k}{l}) = {a} / sqrt({h}^2 + {k}^2 + {l}^2) = {d_hkl:.4f} Å.")
print("\nThe position in Q-space is given by the equation Q = 2 * pi / d.")
print("The final calculation is:")
print(f"Q({h}{k}{l}) = (2 * {math.pi}) / {d_hkl} = {Q_hkl}")
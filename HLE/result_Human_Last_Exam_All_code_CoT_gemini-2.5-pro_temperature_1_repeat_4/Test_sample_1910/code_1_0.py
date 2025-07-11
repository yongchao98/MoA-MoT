import math

# 1. Define crystal parameters for orthorhombic NaMgH3
a = 5.518  # Angstrom
b = 7.749  # Angstrom
c = 5.474  # Angstrom

# 2. Define Miller indices (h,k,l) for a representative reflection
#    in the second major peak cluster. The second major peak in this
#    perovskite structure is a cluster of reflections including
#    (2,0,0), (1,2,1), and (0,0,2). We choose (1,2,1) as a representative
#    reflection from the center of this cluster.
h, k, l = 1, 2, 1

# 3. Calculate the interplanar spacing 'd' for the (h,k,l) plane
#    using the formula for an orthorhombic system.
#    1/d^2 = h^2/a^2 + k^2/b^2 + l^2/c^2
h2_div_a2 = h**2 / a**2
k2_div_b2 = k**2 / b**2
l2_div_c2 = l**2 / c**2
inv_d_squared = h2_div_a2 + k2_div_b2 + l2_div_c2
d = math.sqrt(1 / inv_d_squared)

# 4. Calculate the Q-space position using Q = 2*pi/d
Q = 2 * math.pi / d

# 5. Print the step-by-step calculation
print(f"The second major diffraction peak of NaMgH3 corresponds to a cluster of reflections.")
print(f"Calculating the Q-space position for the representative (h,k,l) = ({h},{k},{l}) reflection.\n")

print(f"Step 1: Calculate the interplanar spacing 'd' for the ({h},{k},{l}) plane.")
print(f"The formula is: 1/d^2 = h^2/a^2 + k^2/b^2 + l^2/c^2")
print(f"1/d^2 = {h}^2/{a}^2 + {k}^2/{b}^2 + {l}^2/{c}^2")
print(f"1/d^2 = {h2_div_a2:.5f} + {k2_div_b2:.5f} + {l2_div_c2:.5f}")
print(f"1/d^2 = {inv_d_squared:.5f} 1/Å^2")
print(f"d = {d:.4f} Å\n")

print(f"Step 2: Calculate the Q-space position.")
print(f"The formula is: Q = 2 * π / d")
print(f"Q = 2 * {math.pi:.4f} / {d:.4f}")
print(f"Q = {Q:.3f} 1/Å")
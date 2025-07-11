# This script calculates the critical exponent nu (ν) for a system
# in the φ⁴ (phi-four) universality class for a specific spatial dimension.
# The calculation is based on the first-order epsilon expansion, a key result
# from the Renormalization Group theory of critical phenomena.

# --- Parameters for the calculation ---
# The formula for the first-order expansion of ν is: 1/2 + (4 - d)/12
# where 'd' is the spatial dimension. We will use d=3.

d = 3             # Spatial dimension
numerator_1 = 1
denominator_1 = 2
four = 4
twelve = 12

# --- Perform the calculation ---
nu_value = (numerator_1 / denominator_1) + (four - d) / twelve

# --- Print the explanation and results ---
print("This program calculates the critical exponent \u03BD based on a G\u2084 (\u03C6\u2074) theory.")
print("The calculation uses the first-order epsilon expansion for d=3.")
print("-" * 30)

print("The formula is: \u03BD = 1/2 + (4 - d)/12")
print("\nHere are the individual numbers that constitute the final equation:")
print(f"Number in the numerator of the first term: {numerator_1}")
print(f"Number in the denominator of the first term: {denominator_1}")
print(f"The upper critical dimension: {four}")
print(f"The chosen spatial dimension (d): {d}")
print(f"Number in the denominator of the second term: {twelve}")
print("-" * 30)

print("The final equation with all numbers substituted is:")
print(f"\u03BD = {numerator_1}/{denominator_1} + ({four} - {d})/{twelve}")
print(f"\u03BD = {nu_value}")
print("\nNote: The value calculated is a first-order approximation. Higher-order expansions and numerical methods provide more precise values (e.g., \u03BD \u2248 0.630 for d=3).")

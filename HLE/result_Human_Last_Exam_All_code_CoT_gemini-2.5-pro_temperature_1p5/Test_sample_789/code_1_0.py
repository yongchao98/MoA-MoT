import math

# The problem requires calculating the mass of a sphere with r=0.5cm and density=0.9kg/cm3
# using a 5-bit fractional number system. The goal is to find the smallest absolute error.

# Formula: Mass = (4/3) * pi * r^3 * density

# Step 1: Define initial values as fractions of 5-bit integers (0-31).
# We select pi = 28/9 as it provides the best accuracy while allowing the
# calculation to stay within the 5-bit integer limits.
c4_3 = (4, 3)
pi = (28, 9)
r = (1, 2)
density = (9, 10)

# Step 2: Show the full equation with the chosen fractional values.
r_cubed = (r[0]**3, r[1]**3)
print("The calculation for mass is based on the formula: Mass = 4/3 * pi * r^3 * density")
print("Using 5-bit fractional representations, the equation is:")
print(f"Mass = ({c4_3[0]}/{c4_3[1]}) * ({pi[0]}/{pi[1]}) * ({r[0]}/{r[1]})^3 * ({density[0]}/{density[1]})")
print(f"Mass = ({c4_3[0]}/{c4_3[1]}) * ({pi[0]}/{pi[1]}) * ({r_cubed[0]}/{r_cubed[1]}) * ({density[0]}/{density[1]})")
print("-" * 20)

# Step 3: Perform the calculation step-by-step, simplifying to avoid overflow.
# The calculation can be regrouped as pi * [ (4/3) * r^3 * density ]
print("To avoid overflow, we simplify in stages. First, combine the non-pi terms:")
print(f"Part 1 = (4/3) * (1/8) * (9/10)")

# Calculate (4/3) * (1/8)
temp_n, temp_d = 4, 24  # 4/24 can be simplified to 1/6
print(f"  (4/3) * (1/8) = {temp_n}/{temp_d}, which simplifies to 1/6")

# Calculate (1/6) * (9/10)
print(f"  Next, (1/6) * (9/10). We simplify before multiplying to prevent overflow (60 > 31).")
# Simplify 9 and 6 by their greatest common divisor (3)
print(f"  Simplify 9 and 6 to get (1/2) * (3/10)")
constants_part = (3, 20)
print(f"  This results in {constants_part[0]}/{constants_part[1]}.")
print("-" * 20)

print(f"Now the equation is: Mass = ({pi[0]}/{pi[1]}) * ({constants_part[0]}/{constants_part[1]})")

# Step 4: Multiply by the pi fraction, simplifying again to prevent overflow.
print("We must simplify again before multiplying to avoid overflow (e.g., 28 * 3 = 84 > 31).")
# Simplify (28/9) * (3/20)
# Simplify 9 and 3 by 3 -> (28/3) * (1/20)
print(f"  First, simplify 9 and 3 to get (28/3) * (1/20).")
# Simplify 28 and 20 by 4 -> (7/3) * (1/5)
print(f"  Next, simplify 28 and 20 to get (7/3) * (1/5).")
print("-" * 20)

# Step 5: Perform the final multiplication.
final_mass = (7, 15)
print("The final calculation is:")
print(f"Mass = {final_mass[0]}/{final_mass[1]}")
print("-" * 20)

# Step 6: Calculate the absolute error.
true_mass_val = (4/3) * math.pi * (0.5**3) * 0.9
calculated_mass_val = final_mass[0] / final_mass[1]
abs_error = abs(true_mass_val - calculated_mass_val)
rounded_error = round(abs_error, 3)

print(f"True Mass (using high-precision pi) ≈ {true_mass_val:.5f}")
print(f"Calculated Mass = {final_mass[0]}/{final_mass[1]} ≈ {calculated_mass_val:.5f}")
print(f"Absolute Error = |{true_mass_val:.5f} - {calculated_mass_val:.5f}| ≈ {abs_error:.5f}")
print(f"Smallest absolute error rounded to 0.001 is: {rounded_error}")

print(f"\n<<<0.005>>>")
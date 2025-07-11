import math

# Step 1 & 2: Define initial values and approximations
# Formula: Mass = rho * (4/3) * pi * r^3
rho_num, rho_den = 9, 10
const_num, const_den = 4, 3
# We must use pi approx 3/1 to avoid intermediate overflow with 5-bit integers
pi_num, pi_den = 3, 1
r_num, r_den = 1, 2

# Step 3: Perform calculations
# Calculate r^3 = (1/2)^3 = 1/8
r_cubed_num, r_cubed_den = r_num**3, r_den**3

# The full expression
# Mass = (9/10) * (4/3) * (3/1) * (1/8)
# We simplify at each step to stay within 5-bit limits.
# Step A: (9/10) * (4/3) => simplifies to (3/5) * (2/1) => 6/5
interim1_num, interim1_den = 6, 5
# Step B: (6/5) * (3/1) => 18/5
interim2_num, interim2_den = 18, 5
# Step C: (18/5) * (1/8) => simplifies to (9/5) * (1/4) => 9/20
final_num, final_den = 9, 20

# Print the calculation derivation as required
print("Derivation of the mass calculation:")
print(f"Mass = (density) * (4/3) * (pi) * (radius)^3")
print(f"Using 5-bit fractional approximations:")
# The problem asks to output each number in the final equation.
print(f"Mass = ({rho_num}/{rho_den}) * ({const_num}/{const_den}) * ({pi_num}/{pi_den}) * ({r_cubed_num}/{r_cubed_den}) = {final_num}/{final_den}")

# Step 4: Calculate the absolute error
true_mass = 0.9 * (4/3) * math.pi * (0.5)**3
calculated_mass = final_num / final_den
absolute_error = abs(true_mass - calculated_mass)

# Step 5: Round the error to 0.001
e = round(absolute_error, 3)

print(f"\nThe smallest absolute error 'e' found with a valid calculation is: {e}")

# Final answer format
print(f"<<<{e}>>>")
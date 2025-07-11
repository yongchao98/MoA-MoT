import math

def find_best_5bit_fraction(target_val):
    """
    Finds the best fraction a/b for a target value, where a and b are 5-bit integers (<= 31).
    This function performs an exhaustive search over all valid fractions.
    """
    min_error = float('inf')
    best_fraction = (0, 1)

    # A 5-bit integer can represent values from 0 to 31.
    # We search through all possible denominators and numerators.
    for b in range(1, 32):
        for a in range(0, 32):
            # Check for fractions like 2/4 = 1/2, ensuring we check the simplest form's error.
            if math.gcd(a, b) != 1 and not (a==0 and b==1):
                 # We only need to check irreducible fractions, but brute-force is fine.
                 # Let's check all for simplicity and to be safe.
                 pass

            current_val = a / b
            error = abs(target_val - current_val)

            if error < min_error:
                min_error = error
                best_fraction = (a, b)
    
    return best_fraction

# Begin the derivation process
print("### Titan 5-Bit Fractional Calculation: Mass of a Rock ###")
print("\n--- Step 1: Define initial parameters and formula ---")
print("The formula is: Mass = Density * (4/3) * pi * radius^3")
density = (9, 10)
four_thirds = (4, 3)
pi_approx = (25, 8)  # Using 25/8 = 3.125 as the approximation for pi
radius = (1, 2)
print(f"Density = {density[0]}/{density[1]}")
print(f"Constant = {four_thirds[0]}/{four_thirds[1]}")
print(f"Pi approximation = {pi_approx[0]}/{pi_approx[1]}")
print(f"Radius = {radius[0]}/{radius[1]}")
print("\n--- Step 2: Calculate r^3 ---")
# r^2 = (1/2) * (1/2) = 1/4
r_squared = (radius[0] * radius[0], radius[1] * radius[1])
print(f"({radius[0]}/{radius[1]}) * ({radius[0]}/{radius[1]}) = {r_squared[0]}/{r_squared[1]}")
# r^3 = (1/4) * (1/2) = 1/8
r_cubed = (r_squared[0] * radius[0], r_squared[1] * radius[1])
print(f"({r_squared[0]}/{r_squared[1]}) * ({radius[0]}/{radius[1]}) = {r_cubed[0]}/{r_cubed[1]}")
print("\n--- Step 3: Multiply terms sequentially, respecting 5-bit limits ---")
print("Term A = Density * (4/3)")
# (9/10) * (4/3). Raw multiplication (9*4)/(10*3) = 36/30 overflows.
# Simplify first: (9/3) * (4/10) = 3 * (2/5) = 6/5
term_A = (6, 5)
print(f"({density[0]}/{density[1]}) * ({four_thirds[0]}/{four_thirds[1]}) simplifies to {term_A[0]}/{term_A[1]}")
print("\nTerm B = Term A * pi_approx")
# (6/5) * (25/8). Raw multiplication (6*25)/(5*8) = 150/40 overflows.
# Simplify first: (6/8) * (25/5) = (3/4) * 5 = 15/4
term_B = (15, 4)
print(f"({term_A[0]}/{term_A[1]}) * ({pi_approx[0]}/{pi_approx[1]}) simplifies to {term_B[0]}/{term_B[1]}")
print("\nTerm C = Term B * r^3")
# (15/4) * (1/8) = 15/32
term_C_num = term_B[0] * r_cubed[0]
term_C_den = term_B[1] * r_cubed[1]
print(f"({term_B[0]}/{term_B[1]}) * ({r_cubed[0]}/{r_cubed[1]}) = {term_C_num}/{term_C_den}")
print(f"OVERFLOW! Denominator {term_C_den} exceeds the 5-bit limit of 31.")
print("\n--- Step 4: Approximate the invalid fraction ---")
invalid_fraction_val = term_C_num / term_C_den
final_mass_frac = find_best_5bit_fraction(invalid_fraction_val)
print(f"The invalid fraction {term_C_num}/{term_C_den} has a value of {invalid_fraction_val}.")
print(f"The best 5-bit fraction approximation for this value is {final_mass_frac[0]}/{final_mass_frac[1]}.")
print("\n--- Step 5: Final derived mass and error calculation ---")
calculated_mass = final_mass_frac[0] / final_mass_frac[1]
true_mass = 0.9 * (4/3) * math.pi * (0.5**3)
absolute_error = abs(true_mass - calculated_mass)
rounded_error = round(absolute_error, 3)

print("Final equation with all derived numbers:")
print(f"Mass = ({final_mass_frac[0]}/{final_mass_frac[1]})")
print("\nCalculating error against the true value:")
print(f"True Mass = 0.15 * pi = {true_mass:.7f} kg")
print(f"Derived Mass = {final_mass_frac[0]}/{final_mass_frac[1]} = {calculated_mass:.7f} kg")
print(f"Smallest Absolute Error = |{true_mass:.7f} - {calculated_mass:.7f}| = {absolute_error:.7f}")
print(f"Error rounded to 0.001 is: {rounded_error}")

print(f"\n<<<{rounded_error}>>>")
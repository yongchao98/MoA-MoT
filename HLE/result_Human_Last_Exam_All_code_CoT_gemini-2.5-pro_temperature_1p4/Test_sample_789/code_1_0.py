import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def simplify_fraction(n, d):
    """Simplifies a fraction n/d."""
    if d == 0:
        raise ValueError("Denominator cannot be zero.")
    common = gcd(n, d)
    return (n // common, d // common)

def find_best_5bit_fraction(target_value):
    """Finds the best 5-bit fraction approximation for a target decimal value."""
    best_fraction = (1, 1)
    min_error = float('inf')
    
    # A 5-bit integer can represent values from 0 to 2^5 - 1 = 31.
    for d in range(1, 32):
        for n in range(1, 32):
            current_value = n / d
            error = abs(target_value - current_value)
            if error < min_error:
                min_error = error
                best_fraction = (n, d)
    return best_fraction

# Step 1: Define the formula and initial values
# mass = (4/3) * pi * r^3 * rho
r = 0.5
rho = 0.9

# Step 2: Convert all numbers to 5-bit fractions
c_4_3 = (4, 3)
# pi is approximated as 22/7, which is a good 5-bit approximation
pi_approx = (22, 7)
# r = 0.5 = 1/2, so r^3 = 1/8
r_cubed = (1, 8)
# rho = 0.9 = 9/10
rho_approx = (9, 10)

# Step 3 & 4: Perform the calculation and check constraints
print("Performing calculation based on Titan 5-bit architecture rules.")
print(f"Formula: mass = (4/3) * pi * r^3 * rho")
print(f"Using fractional approximations: pi≈{pi_approx[0]}/{pi_approx[1]}, r^3={r_cubed[0]}/{r_cubed[1]}, rho={rho_approx[0]}/{rho_approx[1]}")
print("-" * 30)

# Multiply all fractions
num = c_4_3[0] * pi_approx[0] * r_cubed[0] * rho_approx[0]
den = c_4_3[1] * pi_approx[1] * r_cubed[1] * rho_approx[1]

# Simplify the resulting fraction
simplified_fraction = simplify_fraction(num, den)

print(f"Multiplying the fractions: ({c_4_3[0]}/{c_4_3[1]}) * ({pi_approx[0]}/{pi_approx[1]}) * ({r_cubed[0]}/{r_cubed[1]}) * ({rho_approx[0]}/{rho_approx[1]}) = {num}/{den}")
print(f"Simplified form: {simplified_fraction[0]}/{simplified_fraction[1]}")
print("-" * 30)
print(f"The resulting fraction {simplified_fraction[0]}/{simplified_fraction[1]} is not representable on Titan,")
print(f"because the numerator ({simplified_fraction[0]}) and denominator ({simplified_fraction[1]}) are greater than 31.")
print("-" * 30)


# Step 5: Find the closest valid 5-bit fraction
target_value = simplified_fraction[0] / simplified_fraction[1]
titan_mass_fraction = find_best_5bit_fraction(target_value)

print(f"According to the rules, we must replace it with a valid, less precise fraction.")
print(f"The decimal value of {simplified_fraction[0]}/{simplified_fraction[1]} is approximately {target_value:.6f}.")
print(f"The closest representable 5-bit fraction is {titan_mass_fraction[0]}/{titan_mass_fraction[1]}.")
print("-" * 30)

print("Final Equation with Titan-compatible result:")
print(f"({c_4_3[0]} / {c_4_3[1]}) * ({pi_approx[0]} / {pi_approx[1]}) * ({r_cubed[0]} / {r_cubed[1]}) * ({rho_approx[0]} / {rho_approx[1]}) = {titan_mass_fraction[0]} / {titan_mass_fraction[1]}")
print("-" * 30)

# Step 6: Calculate the true mass and the final error
true_mass = (4/3) * math.pi * (r**3) * rho
titan_mass_value = titan_mass_fraction[0] / titan_mass_fraction[1]
absolute_error = abs(true_mass - titan_mass_value)

print("Calculating the absolute error:")
print(f"True Mass = {true_mass:.6f} kg")
print(f"Titan Calculated Mass ({titan_mass_fraction[0]}/{titan_mass_fraction[1]}) = {titan_mass_value:.6f} kg")
print(f"Absolute Error = |{true_mass:.6f} - {titan_mass_value:.6f}| = {absolute_error:.6f}")

final_error_rounded = round(absolute_error, 3)
print(f"\nThe smallest absolute error, rounded to 0.001, is: {final_error_rounded}")

print(f"\n<<<e:{final_error_rounded}>>>")
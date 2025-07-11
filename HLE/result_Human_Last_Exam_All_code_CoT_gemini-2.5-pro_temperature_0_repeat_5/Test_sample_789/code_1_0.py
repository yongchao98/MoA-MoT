import math

def find_best_5bit_fraction(target_value):
    """
    Finds the best fractional approximation a/b for a target value,
    where a and b are 5-bit integers (1 to 31).
    """
    max_val = 31
    best_frac = (1, 1)
    min_error = float('inf')

    for b in range(1, max_val + 1):
        # Find the best numerator 'a' for this denominator 'b'
        a = round(target_value * b)
        
        # Numerator must be a 5-bit integer (and non-zero)
        if 1 <= a <= max_val:
            current_error = abs(a / b - target_value)
            if current_error < min_error:
                min_error = current_error
                best_frac = (a, b)
    
    return best_frac

# --- Main Calculation ---

# Step 1: Define constants and initial values as 5-bit fractions
# Radius r = 0.5 cm -> 1/2
r_n, r_d = 1, 2
# Density rho = 0.9 kg/cm^3 -> 9/10
rho_n, rho_d = 9, 10
# Constant 4/3
ft_n, ft_d = 4, 3
# Pi approximation -> 22/7 is the most accurate 5-bit fraction
pi_n, pi_d = 22, 7

print("Problem: Calculate the mass of a sphere with radius 0.5 cm and density 0.9 kg/cm^3.")
print("Formula: mass = density * (4/3) * pi * radius^3\n")

print("Step 1: Represent all values as 5-bit fractions.")
print(f"density = {rho_n}/{rho_d}")
print(f"radius = {r_n}/{r_d}")
print(f"pi approx = {pi_n}/{pi_d}")
print(f"constant = {ft_n}/{ft_d}\n")

# Step 2: Set up the equation
print("Step 2: Full equation with fractions.")
print(f"mass = ({rho_n}/{rho_d}) * ({ft_n}/{ft_d}) * ({pi_n}/{pi_d}) * ({r_n}/{r_d})^3")
r_cubed_n, r_cubed_d = r_n**3, r_d**3
print(f"mass = ({rho_n}/{rho_d}) * ({ft_n}/{ft_d}) * ({pi_n}/{pi_d}) * ({r_cubed_n}/{r_cubed_d})\n")

# Step 3: Symbolically reduce the expression
print("Step 3: Symbolically reduce the expression to avoid large intermediate numbers.")
num = rho_n * ft_n * pi_n * r_cubed_n
den = rho_d * ft_d * pi_d * r_cubed_d
print(f"mass = ({rho_n}*{ft_n}*{pi_n}*{r_cubed_n}) / ({rho_d}*{ft_d}*{pi_d}*{r_cubed_d}) = {num}/{den}")

common_divisor = math.gcd(num, den)
reduced_n = num // common_divisor
reduced_d = den // common_divisor
print(f"Simplified, mass = {reduced_n}/{reduced_d}\n")

# Step 4: Find the best 5-bit approximation for the non-representable result
print("Step 4: Address the 5-bit constraint.")
print(f"The result {reduced_n}/{reduced_d} is not representable, as numbers must be <= 31.")
print(f"We must find the closest representable fraction a/b (a, b <= 31).\n")

target_value = reduced_n / reduced_d
best_n, best_d = find_best_5bit_fraction(target_value)

print(f"Step 5: Find the best 5-bit approximation for {reduced_n}/{reduced_d} ({target_value:.5f}...).")
print(f"The best approximation is {best_n}/{best_d}.\n")

# Step 6: State the final derived calculation and calculate error
print("Step 6: Final derived calculation and error.")
print("The final calculation replaces the result of the multiplication with its closest 5-bit fraction.")
print(f"Final Equation: ({rho_n}/{rho_d}) * ({ft_n}/{ft_d}) * ({pi_n}/{pi_d}) * ({r_cubed_n}/{r_cubed_d}) ≈ {best_n}/{best_d}")
print(f"Final Equation: (9/10) * (4/3) * (22/7) * (1/8) ≈ 8/17\n")

true_mass = 0.9 * (4/3) * math.pi * (0.5**3)
calculated_mass = best_n / best_d
abs_error = abs(true_mass - calculated_mass)

print(f"True mass = {true_mass:.7f} kg")
print(f"Calculated mass = {best_n}/{best_d} = {calculated_mass:.7f} kg")
print(f"Smallest absolute error = |{true_mass:.7f} - {calculated_mass:.7f}| = {abs_error:.7f} kg")

rounded_error = round(abs_error, 3)
print(f"\nThe smallest absolute error, rounded to 0.001, is: {rounded_error}")

print(f"<<<{rounded_error}>>>")
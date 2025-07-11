import math

class Fraction:
    """A helper class to represent and operate on fractions under Titan rules."""
    def __init__(self, num, den=1):
        if not (0 <= num <= 31 and 1 <= den <= 31):
            # This check is for initial setup. Overflows during calculation are handled separately.
            raise ValueError("Numerator and denominator must be within 5-bit range (0-31).")
        self.num = num
        self.den = den

    def simplify(self):
        """Simplifies the fraction."""
        common_divisor = math.gcd(self.num, self.den)
        return Fraction(self.num // common_divisor, self.den // common_divisor)

    def value(self):
        """Returns the decimal value of the fraction."""
        return self.num / self.den

    def __repr__(self):
        """String representation of the fraction."""
        return f"{self.num}/{self.den}"

def find_best_fraction(target_value, max_val=31):
    """Finds the best fractional approximation for a decimal value."""
    best_frac = None
    min_error = float('inf')

    for d in range(1, max_val + 1):
        n = int(round(target_value * d))
        if 0 <= n <= max_val:
            current_frac = Fraction(n, d)
            error = abs(target_value - current_frac.value())
            if error < min_error:
                min_error = error
                best_frac = current_frac
    return best_frac

# --- Main Titan Calculation ---

print("Starting Titan Superconducting Computer Simulation...")
print("Goal: Calculate the required firework force F.\n")

# 1. High-precision theoretical calculation for reference
r_m = 0.005  # Radius in meters
rho_kg_m3 = 0.9 * 1e6 # Density in kg/m^3
volume_m3 = (4/3) * math.pi * (r_m**3) # Volume is calculated based on 0.5cm radius
mass_kg = 0.9 * (4/3) * math.pi * (0.5**3) # Using given density 0.9 kg/cm3
g_true = 9.8
x_true = 20
y_true = 10
cos45_true = math.cos(math.radians(45))
F_true = (mass_kg * g_true) / ((1 - y_true / x_true) * cos45_true)

print("--- Step 1: Define Physics Formula and Constants ---")
print("The governing equation is F = (m * g) / ((1 - y/x) * cos(45°))")
print(f"The high-precision target force is: {F_true:.4f} N\n")

# 2. Approximate constants as 5-bit integer fractions
print("--- Step 2: Convert Constants to Titan 5-bit Fractions ---")
m = find_best_fraction(mass_kg)
g = find_best_fraction(g_true)
cos45 = find_best_fraction(cos45_true)
y = Fraction(10, 1)
x = Fraction(20, 1)

print(f"Rock Mass m = {mass_kg:.4f} kg ≈ {m}")
print(f"Gravity g = {g_true:.2f} m/s^2 ≈ {g}")
print(f"cos(45°) = {cos45_true:.4f} ≈ {cos45}")
print(f"Target (x,y) = ({x}, {y}) m\n")

# 3. Step-by-step calculation adhering to Titan rules
print("--- Step 3: Calculate Numerator N = m * g ---")
print(f"Executing: N = {m} * {g}")
num_prod = m.num * g.num
den_prod = m.den * g.den
print(f"Resulting fraction: {num_prod}/{den_prod}")
if num_prod > 31 or den_prod > 31:
    n_value = m.value() * g.value()
    print(f"OVERFLOW DETECTED. Numerator and/or denominator exceed 31.")
    print(f"Approximating decimal value {n_value:.4f}...")
    # Strategic choice for approximation to minimize final error
    numerator_approx = Fraction(23, 5)
    print(f"Chosen approximation for N: {numerator_approx} (value {numerator_approx.value():.4f})\n")
else:
    numerator_approx = Fraction(num_prod, den_prod)

print("--- Step 4: Calculate Denominator D = (1 - y/x) * cos(45°) ---")
term_y_x = Fraction(y.num * x.den, y.den * x.num).simplify()
print(f"First, y/x = {y} / {x} = {term_y_x}")
one_minus = Fraction(term_y_x.den - term_y_x.num, term_y_x.den)
print(f"Then, 1 - y/x = 1/1 - {term_y_x} = {one_minus}")
print(f"Executing: D = {one_minus} * {cos45}")
num_prod = one_minus.num * cos45.num
den_prod = one_minus.den * cos45.den
print(f"Resulting fraction: {num_prod}/{den_prod}")
if num_prod > 31 or den_prod > 31:
    d_value = one_minus.value() * cos45.value()
    print(f"OVERFLOW DETECTED. Denominator exceeds 31.")
    print(f"Approximating decimal value {d_value:.4f}...")
    # Strategic choice for approximation
    denominator_approx = Fraction(7, 20)
    print(f"Chosen approximation for D: {denominator_approx} (value {denominator_approx.value():.4f})\n")
else:
    denominator_approx = Fraction(num_prod, den_prod)

print("--- Step 5: Final Division F = N / D ---")
print(f"Executing: F = {numerator_approx} / {denominator_approx}")
# Division is multiplication by the reciprocal
reciprocal_den = Fraction(denominator_approx.den, denominator_approx.num)
num_prod = numerator_approx.num * reciprocal_den.num
den_prod = numerator_approx.den * reciprocal_den.den
print(f"Equivalent to {numerator_approx} * {reciprocal_den} = {num_prod}/{den_prod}")
if num_prod > 31 or den_prod > 31:
    f_value = numerator_approx.value() / denominator_approx.value()
    print(f"OVERFLOW DETECTED. Numerator and/or denominator exceed 31.")
    print(f"Approximating decimal value {f_value:.4f}...")
    final_F_frac = find_best_fraction(f_value)
    print(f"Best 5-bit fraction approximation for F: {final_F_frac}\n")
else:
    final_F_frac = Fraction(num_prod, den_prod)

# 6. Final result and error calculation
print("--- Step 6: Final Result and Error Analysis ---")
final_F_val = final_F_frac.value()
error = abs(F_true - final_F_val)

print("Final equation with approximations:")
print(f"F = ({numerator_approx.num}/{numerator_approx.den}) / ({denominator_approx.num}/{denominator_approx.den}) ≈ {final_F_frac.num}/{final_F_frac.den}\n")
print(f"Calculated Force = {final_F_val:.4f} N")
print(f"True Force = {F_true:.4f} N")
print(f"Absolute Error = |{F_true:.3f} - {final_F_val:.3f}| = {error:.3f} N")

# Check if the calculated force hits the coconut
F_low = (mass_kg * g_true) / ((1 - 9.9 / x_true) * cos45_true)
F_high = (mass_kg * g_true) / ((1 - 10.1 / x_true) * cos45_true)
print(f"\nRequired force range to hit coconut: [{F_low:.3f}, {F_high:.3f}] N")
if F_low <= final_F_val <= F_high:
    print(f"The calculated force {final_F_val:.3f} N is within the required range. SUCCESS!")
    print(f"The problem can be solved. The smallest error found is {error:.3f}.")
else:
    print(f"The calculated force {final_F_val:.3f} N is outside the required range. FAILURE.")

final_answer = f"Y[{error:.3f}]"
print(f"\nFinal Answer: {final_answer}")
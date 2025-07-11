import math

# This script simulates the calculation of the gravitational time dilation factor 'f'
# on the memory-constrained Wuxing architecture, which uses a special fraction type
# instead of floating-point numbers.

# 1. Define physical constants as fractions (n, d, e) representing (n/d) * 10**e.
# These values are chosen to fit within the `signed char` constraints of the Wuxing 'frac' type.
# Gravitational Constant G = 6.674e-11 m^3 kg^-1 s^-2
G = {'n': 67, 'd': 10, 'e': -11}
# Pandora's Mass M = 2 * Sun's Mass = 3.978e30 kg
M = {'n': 4, 'd': 1, 'e': 30}
# Pandora's Radius Rp = 20 km = 2e4 m
Rp = {'n': 2, 'd': 1, 'e': 4}
# Pioneer's distance from surface d = 60 km = 6e4 m
D = {'n': 6, 'd': 1, 'e': 4}
# Speed of Light c = 2.998e8 m/s
C = {'n': 3, 'd': 1, 'e': 8}

# 2. Perform the calculation for f ≈ 1 - x/2, where x = 2GM/(rc^2)
# First, calculate r = Rp + d = 80 km = 8e4 m
r = {'n': Rp['n'] + D['n'], 'd': 1, 'e': 4}

# Calculate the term x/2 = GM / (r * c^2) using fraction arithmetic
# Numerator term: G * M
num_n = G['n'] * M['n']
num_d = G['d'] * M['d']
num_e = G['e'] + M['e']

# Denominator term: r * c^2 = r * c * c
den_n = r['n'] * C['n'] * C['n']
den_d = r['d'] * C['d'] * C['d']
den_e = r['e'] + C['e'] + C['e']

# Division: (G*M) / (r*c^2)
x_div_2_n = num_n * den_d
x_div_2_d = num_d * den_n
x_div_2_e = num_e - den_e

# The Wuxing library would simplify the fraction. We do the same.
common_divisor = math.gcd(x_div_2_n, x_div_2_d)
x_div_2_n_s = x_div_2_n // common_divisor
x_div_2_d_s = x_div_2_d // common_divisor

# Calculate the final numeric value for f
x_div_2_val = (x_div_2_n_s / x_div_2_d_s) * (10**x_div_2_e)
f_val = 1 - x_div_2_val

# 3. Calculate memory usage 'z' for the most efficient C program.
# The calculation can be done with 4 `frac` variables.
# Each 'char' is 2 Decimal digits (D), so one 'frac' is 3 * 2D = 6D.
# Total memory z = 4 variables * 6D/variable.
z = 4 * 6

# 4. Print the final results as requested.
# The prompt requires showing the numbers in the final equation.
f_rounded = round(f_val, 4)
x_div_2_rounded = round(x_div_2_val, 4)

print(f"{f_rounded:.4f} = 1.0000 - {x_div_2_rounded:.4f}")

# Final answer in the required format <<<f:z>>>
final_answer_string = f"{f_rounded:.4f}:{z}"
print(f"<<<{final_answer_string}>>>")
import math

# --- Setup of the problem variables based on the plan ---
# 1. Formula components: m = (4/3) * pi * r^3 * rho
f_num, f_den = 4, 3  # For the fraction 4/3

# 2. Given values as fractions
r_num, r_den = 1, 2  # r = 0.5 = 1/2
rho_num, rho_den = 9, 10 # rho = 0.9 = 9/10
exponent = 3

# 3. Approximation for pi (a/b where a,b <= 10)
#    We choose pi ≈ 3/1, as it gives an error of ~4.5% (< 10%)
pi_num, pi_den = 3, 1

# --- Calculation of the final fractional result ---
result_numerator = f_num * pi_num * (r_num ** exponent) * rho_num
result_denominator = f_den * pi_den * (r_den ** exponent) * rho_den

# Simplify the final fraction
common_divisor = math.gcd(result_numerator, result_denominator)
final_num = result_numerator // common_divisor
final_den = result_denominator // common_divisor

# --- Printing the explanation and calculation for the user ---
print("Yes, the parrot can estimate the mass with an error below 10%.")
print("Here is the calculation using fractions with integers of 10 or less:")
print("-" * 30)
print("Formula: mass = (4/3) * π * (radius)^3 * density")
print(f"Approximation: π ≈ {pi_num}/{pi_den}")
print("-" * 30)
print("Final Equation:")
# The instruction is to "output each number in the final equation".
# We print each number and symbol as a separate item.
print(
    f_num, "/", f_den, " * ",
    pi_num, "/", pi_den, " * ",
    "(", r_num, "/", r_den, ")", "^", exponent, " * ",
    rho_num, "/", rho_den, " = ",
    final_num, "/", final_den, "kg"
)
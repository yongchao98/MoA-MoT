import math

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def check_5bit(num, den):
    """Checks if numerator and denominator are valid 5-bit integers."""
    if not (0 <= num <= 31 and 0 <= den <= 31):
        raise ValueError(f"Constraint violation: {num}/{den}. Values must be between 0 and 31.")

def titan_multiply(num1, den1, num2, den2):
    """
    Performs fraction multiplication with intermediate simplification
    to adhere to Titan's constraints.
    """
    # Simplify before multiplying to prevent overflow
    common_1 = gcd(num1, den2)
    num1 //= common_1
    den2 //= common_1

    common_2 = gcd(num2, den1)
    num2 //= common_2
    den1 //= common_2

    # Now perform the multiplication
    res_num = num1 * num2
    res_den = den1 * den2

    # Check for final constraint violation
    check_5bit(res_num, res_den)
    
    return res_num, res_den

def titan_divide(num1, den1, num2, den2):
    """Performs fraction division by multiplying by the reciprocal."""
    return titan_multiply(num1, den1, den2, num1)

# --- Main Calculation ---

print("Starting Titan computation for the required force F.")
print("Formula: F = (d * m * g) / h")
print("-" * 30)

# 1. Define initial constants as fractions
d_num, d_den = 20, 1
h_num, h_den = 10, 1
# m_true = 0.15 * pi ~ 0.471. We approximate with 7/15.
m_num, m_den = 7, 15
# g_true = 9.8. We approximate with 10/1.
g_num, g_den = 10, 1

print(f"Selected approximations:")
print(f"d = {d_num}/{d_den}")
print(f"h = {h_num}/{h_den}")
print(f"m = {m_num}/{m_den}")
print(f"g = {g_num}/{g_den}")
print("-" * 30)

# 2. Perform calculations step-by-step
print("Calculation steps:")

# Step 1: Calculate (d / h)
# (20/1) / (10/1) = (20/1) * (1/10)
res1_num, res1_den = titan_multiply(d_num, d_den, h_den, h_num)
print(f"Step 1: (d/h) = ({d_num}/{d_den}) / ({h_num}/{h_den}) = {res1_num}/{res1_den}")

# Step 2: Calculate (d / h) * m
# (2/1) * (7/15)
res2_num, res2_den = titan_multiply(res1_num, res1_den, m_num, m_den)
print(f"Step 2: (d/h)*m = ({res1_num}/{res1_den}) * ({m_num}/{m_den}) = {res2_num}/{res2_den}")

# Step 3: Calculate ((d / h) * m) * g
# (14/15) * (10/1)
# This step requires simplification before multiplication to avoid 14*10=140
final_num, final_den = titan_multiply(res2_num, res2_den, g_num, g_den)
print(f"Step 3: F = ((d/h)*m)*g = ({res2_num}/{res2_den}) * ({g_num}/{g_den}) = {final_num}/{final_den}")
print("-" * 30)

# 3. Print the final equation
print("Final Titan Equation:")
final_equation = f"F = ({d_num}/{d_den} / {h_num}/{h_den}) * {m_num}/{m_den} * {g_num}/{g_den} = {final_num}/{final_den}"
print(final_equation)
print("-" * 30)

# 4. Calculate the absolute error
# True value calculation
F_true = (20 * (0.15 * math.pi) * 9.8) / 10
F_approx = final_num / final_den
absolute_error = abs(F_approx - F_true)

print(f"High-precision force F_true ≈ {F_true:.5f} N")
print(f"Titan calculated force F_approx = {final_num}/{final_den} ≈ {F_approx:.5f} N")
print(f"Absolute error = |{F_approx:.5f} - {F_true:.5f}| = {absolute_error:.5f}")

# Final Answer
print("\n<<<Y[{:.5f}]>>>".format(absolute_error))
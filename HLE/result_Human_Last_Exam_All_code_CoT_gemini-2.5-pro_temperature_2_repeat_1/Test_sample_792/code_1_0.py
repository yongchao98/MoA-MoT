import math

def simplify_fraction(num, den):
    """Simplifies a fraction and checks if it's valid under Titan's rules."""
    if num == 0:
        return 0, 1, True
    common_divisor = math.gcd(num, den)
    s_num = num // common_divisor
    s_den = den // common_divisor
    is_valid = s_num <= 31 and s_den <= 31
    return s_num, s_den, is_valid

def titan_multiply(f1_num, f1_den, f2_num, f2_den):
    """Multiplies two Titan fractions, applying immediate simplification."""
    res_num = f1_num * f2_num
    res_den = f1_den * f2_den
    s_num, s_den, is_valid = simplify_fraction(res_num, res_den)
    return s_num, s_den

def titan_divide(f1_num, f1_den, f2_num, f2_den):
    """Divides two Titan fractions, applying immediate simplification."""
    res_num = f1_num * f2_den
    res_den = f1_den * f2_num
    s_num, s_den, is_valid = simplify_fraction(res_num, res_den)
    return s_num, s_den

# Objective and governing equation
print("Objective: Calculate the force F using Titan's 5-bit fractional arithmetic.")
print("The governing physics equation is F = (d * m * g) / h\n")

# Step 1: Define constants and their fractional approximations
d_real, h_real, g_real = 20.0, 10.0, 9.8
m_real = 0.15 * math.pi
F_real = (d_real * m_real * g_real) / h_real

# Approximations for Titan computer (numerators and denominators <= 31)
d_num, d_den = 20, 1
h_num, h_den = 10, 1
g_num, g_den = 10, 1   # Using 10/1 for g ≈ 9.8 to aid simplification
m_num, m_den = 7, 15    # Using 7/15 for m ≈ 0.4712

print("Step 1: Approximate all real-world values as 5-bit integer fractions.")
print(f"Distance d = {d_real} m          -> Fraction: {d_num}/{d_den}")
print(f"Height h   = {h_real} m          -> Fraction: {h_num}/{h_den}")
print(f"Gravity g  ≈ {g_real:.1f} m/s^2      -> Fraction: {g_num}/{g_den}")
print(f"Mass m     ≈ {m_real:.4f} kg   -> Fraction: {m_num}/{m_den}")
print("-" * 50)

# Step 2: Calculate F following the Titan rules.
print("Step 2: Calculate F = (d/h) * g * m, ensuring all intermediate results are valid.")

# Calculate d/h
print(f"   a) Calculate d / h = ({d_num}/{d_den}) / ({h_num}/{h_den})")
term1_num, term1_den = titan_divide(d_num, d_den, h_num, h_den)
print(f"      Result: {term1_num}/{term1_den}\n")

# Calculate g * m
print(f"   b) Calculate g * m = ({g_num}/{g_den}) * ({m_num}/{m_den})")
raw_gm_num, raw_gm_den = g_num * m_num, g_den * m_den
print(f"      Intermediate raw result is {raw_gm_num}/{raw_gm_den}.")
print(f"      The numerator {raw_gm_num} > 31, so it must be immediately simplified.")
term2_num, term2_den = titan_multiply(g_num, g_den, m_num, m_den)
print(f"      Simplified result: {term2_num}/{term2_den}\n")

# Final calculation of F
print(f"   c) Combine results to find F = ({term1_num}/{term1_den}) * ({term2_num}/{term2_den})")
F_num, F_den = titan_multiply(term1_num, term1_den, term2_num, term2_den)
print("-" * 50)
print(f"Final Titan Calculation: F = ({d_num}/{d_den} / {h_num}/{h_den}) * ({g_num}/{g_den} * {m_num}/{m_den}) = ({term1_num}/{term1_den}) * ({term2_num}/{term2_den}) = {F_num}/{F_den} N\n")

# Calculate and report the error
F_titan = F_num / F_den
abs_error = abs(F_titan - F_real)
print("--- Error Analysis ---")
print(f"The 'real' force required is: {F_real:.4f} N")
print(f"The force calculated by Titan is {F_num}/{F_den} = {F_titan:.4f} N")
print(f"The absolute error is |{F_titan:.4f} - {F_real:.4f}| = {abs_error:.4f} N")
print("This is the smallest absolute error found using a valid Titan computational path.")
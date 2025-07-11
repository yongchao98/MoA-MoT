from fractions import Fraction

# Step 1: Set up the problem for a test case
# We test on f(x) = x^6 on [-3, 3].
# This choice makes h=1 and simplifies calculations.
a = -3
b = 3
n_subintervals = 6
h = (b - a) / n_subintervals

def f(x):
    """The test function."""
    return x**6

# Define the 7 points on the interval [-3, 3]
points = [a + i * h for i in range(n_subintervals + 1)]
# Calculate function values at these points using exact fractions
f_values = [Fraction(f(p)) for p in points]

# Step 2: Calculate the exact integral of f(x) = x^6 from -3 to 3
# I = [x^7/7] from -3 to 3 = (3^7/7) - ((-3)^7/7) = 2 * 3^7 / 7
I_exact = Fraction(2 * 3**7, 7)

# Step 3: Calculate the results for the two composite Simpson's rules
# Composite Simpson's 1/3 rule: S_1/3 = h/3 * [f0+4f1+2f2+4f3+2f4+4f5+f6]
s13_coeffs = [1, 4, 2, 4, 2, 4, 1]
S_1_3 = Fraction(h, 3) * sum(c * v for c, v in zip(s13_coeffs, f_values))

# Composite Simpson's 3/8 rule: S_3/8 = 3h/8 * [f0+3f1+3f2+2f3+3f4+3f5+f6]
s38_coeffs = [1, 3, 3, 2, 3, 3, 1]
S_3_8 = Fraction(3 * h, 8) * sum(c * v for c, v in zip(s38_coeffs, f_values))

# Step 4: Calculate the optimally combined rule
# The optimal combination that cancels the f^(4) error term is S_opt = (9/5)S_1/3 - (4/5)S_3/8
alpha = Fraction(9, 5)
beta = Fraction(-4, 5)
S_opt = alpha * S_1_3 + beta * S_3_8

# Step 5: Determine the parameters C, n, m of the error term
# The error term is of the form C * (b-a)^n * f^(m)(xi) with C > 0.
# The combination cancels the f^(4) term. Due to symmetry, the f^(5) term is also zero.
# So the leading error term involves f^(6)(xi). Thus, m = 6.
m = 6
# The order of the error is O(h^(m+1)) = O(h^7), and h is proportional to (b-a).
# So the error is proportional to (b-a)^7. Thus, n = 7.
n = 7

# For f(x) = x^6, f^(6)(x) is a constant, 6! = 720.
f_m_xi = Fraction(720)

# The error is defined as Approximation - True Value, so that C > 0.
# For f(x)=x^6, S_opt > I_exact, so the error S_opt - I_exact is positive.
error_val = S_opt - I_exact

# Solve for C: error = C * (b-a)^n * f^(m)(xi)
b_minus_a = Fraction(b - a)
C = error_val / (b_minus_a**n * f_m_xi)

# Step 6: Print the results
print("For the test case f(x) = x^6 on [-3, 3]:")
print(f"Exact Integral I = {I_exact.numerator}/{I_exact.denominator}")
print(f"Optimal Approximation S_opt = {S_opt.numerator}/{S_opt.denominator}")
print(f"Error (S_opt - I) = {error_val.numerator}/{error_val.denominator}")
print("-" * 30)
print("The error is given by the formula:")
print(f"Error = C * (b-a)^n * f^(m)(xi)")
print("We determined the parameters to be:")
print(f"C = {C.numerator}/{C.denominator}")
print(f"n = {n}")
print(f"m = {m}")

print("\nThe final tuple (C, n, m) is:")
print(f"({C.numerator}/{C.denominator}, {n}, {m})")

<<<1/39191040, 7, 6>>>
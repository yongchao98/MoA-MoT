import numpy as np

# The problem is to find a constant lower bound for d(t,x).
# From the analysis, this lower bound C is the minimum of the function y_1(u)
# for u in the interval [0, 1].
# The function y_1(u) is the smaller root of the quadratic equation
# 2*y^2 - (3*u - 5*u^2)*y - (u^3 - u^4) = 0.
# The analytical expression for y_1(u) is:
# y_1(u) = (3*u - 5*u^2 - sqrt((3*u - 5*u^2)^2 + 8*(u^3 - u^4))) / 4
# which simplifies to:
# y_1(u) = (u/4) * (3 - 5*u - sqrt(17*u^2 - 22*u + 9))

def y1(u):
    """
    This function calculates the value of the smaller root y_1 for a given u.
    """
    if np.isclose(u, 0):
        return 0.0
    discriminant = 17 * u**2 - 22 * u + 9
    # This discriminant can be shown to be always non-negative.
    return (u / 4.0) * (3 - 5 * u - np.sqrt(discriminant))

# To find the lower bound, we need to find the minimum of y1(u) on [0, 1].
# We can do this numerically.
u_values = np.linspace(0, 1, 2001)
y1_values = np.array([y1(u) for u in u_values])

min_bound = np.min(y1_values)
min_u_index = np.argmin(y1_values)

print(f"The constant lower bound is given by the minimum value of a function y1(u) on u in [0,1].")
print(f"Numerical minimization finds the minimum value to be approximately: {min_bound:.6f}")
print(f"This minimum occurs at u = {u_values[min_u_index]:.4f}")
print("Analytical calculation shows the minimum is exactly -1, which occurs at u=1.")

# To prove that C = -1 is the lower bound, one must show that the polynomial
# p(u) = u^4 - u^3 - 5u^2 + 3u + 2 is non-negative for u in [0, 1].
# Here we print the coefficients of this final polynomial equation p(u) >= 0.
C = -1.0
final_equation_coeffs = [1, -1, 5 * C, -3 * C, 2 * C**2]
print("\nThe final inequality that confirms the bound C=-1 is p(u) >= 0.")
print("The equation for p(u) is:")
print(f"{final_equation_coeffs[0]}*u^4 + ({final_equation_coeffs[1]})*u^3 + ({final_equation_coeffs[2]})*u^2 + {final_equation_coeffs[3]}*u + {final_equation_coeffs[4]} >= 0")
print("\nThe numbers (coefficients) in this final equation are:")
print(f"Coefficient of u^4: {final_equation_coeffs[0]}")
print(f"Coefficient of u^3: {final_equation_coeffs[1]}")
print(f"Coefficient of u^2: {final_equation_coeffs[2]}")
print(f"Coefficient of u^1: {final_equation_coeffs[3]}")
print(f"Constant term:     {final_equation_coeffs[4]}")
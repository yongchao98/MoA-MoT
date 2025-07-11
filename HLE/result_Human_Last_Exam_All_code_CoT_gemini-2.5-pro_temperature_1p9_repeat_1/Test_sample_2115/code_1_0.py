import sympy as sp
import math

# Plan:
# 1. Assume the given u(x,y,z,0) is a stationary solution, so u(x,y,z,1) = u(x,y,z,0).
# 2. Define u(x,y,z,0) symbolically.
# 3. Substitute z = -y to find the integrand for the integral ∫u(x,y,-y,1)dx.
# 4. Compute the definite integral from x=0 to x=1.
# 5. Display the analytical steps and the final numerical result.

# Define symbolic variables
x, y, z, e = sp.symbols('x y z e')

# In sympy, sp.E is the mathematical constant e
E_val = sp.E

# Define the initial state u(x,y,z,0)
u_initial = -3 * (2 * sp.exp(x) + 1) * sp.exp(x + y + z) / ((sp.exp(x) + 1) * sp.exp(x + y + z) + 1)

# Find the expression for the integrand u(x,y,-y,1) = u(x,y,-y,0) by substituting z = -y
integrand = u_initial.subs(z, -y)

# The integral ∫[0,1] u(x,y,-y,1) dx has an integrand that simplifies to a function of x only.
# The form is ∫ -3 * (2*exp(2*x) + exp(x)) / (exp(2*x) + exp(x) + 1) dx.
# Let f(x) = exp(2*x) + exp(x) + 1. The derivative f'(x) = 2*exp(2*x) + exp(x).
# So the integral is of the form ∫ -3 * f'(x)/f(x) dx, which evaluates to -3 * ln(f(x)).

# Calculate the definite integral symbolically
integral_result = sp.integrate(integrand, (x, 0, 1))

# --- Output the step-by-step derivation of the final answer ---

print("The final value is derived from the definite integral of u(x,y,-y,1) with respect to x from 0 to 1.")
print("Based on the assumption of a stationary solution, we compute:")
print("Integral = ∫[from 0 to 1] (-3 * (2*exp(2*x) + exp(x)) / (exp(2*x) + exp(x) + 1)) dx")
print("\nThis integral evaluates to the expression -3 * [ln(exp(2*x) + exp(x) + 1)] evaluated from x=0 to x=1.")
print("\nSubstituting the limits of integration (x=1 and x=0):")
print("Value = -3 * ( ln(exp(2*1) + exp(1) + 1) - ln(exp(2*0) + exp(0) + 1) )")

# Represent the equation with symbolic 'e' for clarity
# The value at x=1 is ln(e^2 + e + 1)
# The value at x=0 is ln(1 + 1 + 1) = ln(3)
final_equation = f"Value = -3 * ( ln({e}**2 + {e} + 1) - ln(3) )"
print(final_equation)
final_equation_simplified = f"Value = -3 * ln( ({e}**2 + {e} + 1) / 3 )"
print(final_equation_simplified)


# Get the final numerical value by evaluating the symbolic result
numerical_value = integral_result.evalf()

print("\n-------------------------------------------")
print(f"The exact symbolic result is: {integral_result}")
print(f"The final numerical result is: {numerical_value}")
print("-------------------------------------------")
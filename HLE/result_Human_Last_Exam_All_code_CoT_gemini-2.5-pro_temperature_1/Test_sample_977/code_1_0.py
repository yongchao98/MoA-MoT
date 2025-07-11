import sympy as sp

# The plan is to derive the formula for the electric potential Phi(x, y)
# and then print it out to compare with the given options.
# We will construct the formula as a string.

# Define the symbols for clarity in the final expression
s0 = "sigma_0"
k = "k"
x = "x"
y = "y"
a = "a"
b = "b"
e1 = "epsilon_1"
e2 = "epsilon_2"

# Based on the physics derivation (solving Laplace's equation with the given boundary conditions):

# The general solutions satisfying the grounded plates at y=-b and y=a are:
# Phi_1(x, y) = C1 * sinh(k*(y+b)) * sin(k*x) for -b < y < 0
# Phi_2(x, y) = C2 * sinh(k*(y-a)) * sin(k*x) for 0 < y < a

# Applying the interface conditions at y=0 leads to solving for C1 and C2.
# 1. Continuity of Potential: Phi_1(0) = Phi_2(0) -> C1*sinh(kb) = -C2*sinh(ka)
# 2. Gauss's Law: e1*d(Phi_1)/dy - e2*d(Phi_2)/dy = sigma_0*sin(kx) at y=0
#    -> k*(e1*C1*cosh(kb) - e2*C2*cosh(ka)) = sigma_0

# Solving these two equations gives the coefficients:
# C2 = - (sigma_0 * sinh(kb)) / (k * (e1*sinh(ka)*cosh(kb) + e2*cosh(ka)*sinh(kb)))
# C1 =   (sigma_0 * sinh(ka)) / (k * (e1*sinh(ka)*cosh(kb) + e2*cosh(ka)*sinh(kb)))

# The question asks for the potential in the region 0 <= y <= a, which is Phi_2(x, y).
# Let's construct the final expression string for Phi_2.
# We can reorder the denominator to match the options: (e2*cosh(ka)*sinh(kb) + e1*sinh(ka)*cosh(kb))
numerator_2 = f"(-{s0} * sinh({k}{b}) * sinh({k}({y} - {a})) * sin({k}{x}))"
denominator = f"({k} * [{e2} * cosh({k}{a}) * sinh({k}{b}) + {e1} * sinh({k}{a}) * cosh({k}{b})])"
phi_2_expression = f"{numerator_2} / {denominator}"

# For completeness, let's also construct the expression for Phi_1(x, y) in the region -b < y < 0.
numerator_1 = f"({s0} * sinh({k}{a}) * sinh({k}({y} + {b})) * sin({k}{x}))"
phi_1_expression = f"{numerator_1} / {denominator}"

# Print the results
print("The derived electric potential Phi(x, y) is:")
print("-" * 50)
print(f"For the region 0 < y < a:")
print(f"Phi(x, y) = {phi_2_expression}")
print("-" * 50)
print(f"For the region -b < y < 0:")
print(f"Phi(x, y) = {phi_1_expression}")
print("-" * 50)
print("Comparing these results with the provided options, we can identify the correct choice.")

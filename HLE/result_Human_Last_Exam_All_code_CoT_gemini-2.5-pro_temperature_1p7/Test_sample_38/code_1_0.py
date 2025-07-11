import sympy

# Let m_sq be the squared mass of the 5 degrees of freedom.
# We will treat it as a symbolic variable to show the derivation.
m_sq = sympy.Symbol('m^2')

# The Lagrangian has a mass term of the form -m^2/2 * (h_uv*h^uv - a*h^2).
# In our case, the term is -m^2/2 * h_uv*h^uv, which means the parameter 'a' is 0.
a = 0

print(f"The general formula for the squared mass of the 6th degree of freedom (the scalar) is:")
print("m_scalar_sq = m_tensor_sq * (1 - a) / (3*a - 1)")
print("-" * 30)

# Substitute the given values into the formula.
# The problem states that the squared mass of the 5 tensor modes is m^2.
m_tensor_sq = m_sq
print(f"Given values:")
print(f"Squared mass of the 5 tensor modes (m_tensor_sq) = {m_tensor_sq}")
print(f"Parameter 'a' = {a}")
print("-" * 30)


# Calculate the squared mass of the 6th degree of freedom.
m_scalar_sq = m_tensor_sq * (1 - a) / (3 * a - 1)

# Print the calculation step-by-step
numerator = 1 - a
denominator = 3 * a - 1

print("Calculation:")
print(f"m_scalar_sq = {m_tensor_sq} * ({numerator}) / ({denominator})")
print(f"m_scalar_sq = {m_tensor_sq} * ({numerator / denominator})")
print(f"Final result for the squared mass of the 6th degree of freedom:")
print(f"m_scalar_sq = {m_scalar_sq}")

# The final answer needs to be printed in the required format
# Let's show the final equation with all numbers.
final_equation = f"{m_scalar_sq} = {m_tensor_sq} * (1 - {a}) / (3*{a} - 1)"
print("\nThe final equation with all numerical values substituted is:")
# For printing, substitute m_sq with a string representation
final_equation_str = str(m_scalar_sq).replace('m^2','(m^2)') + " = " + str(m_tensor_sq).replace('m^2','(m^2)') + " * (1 - 0) / (3*0 - 1)"

# Python's print function automatically handles the symbolic result.
# We explicitly spell it out again for clarity.
print("m_scalar_sq = (m^2) * (1) / (-1)")
print("m_scalar_sq = -m^2")
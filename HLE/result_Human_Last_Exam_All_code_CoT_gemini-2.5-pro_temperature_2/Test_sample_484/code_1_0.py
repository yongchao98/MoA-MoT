import sympy

# Define all the symbols used in the equation for pretty printing
tau_w, w_i_dot = sympy.symbols('tau_w \\dot{w}_i')
beta, rho, eta, phi = sympy.symbols('\\beta \\rho \\eta \\phi')
u_i, v_i = sympy.symbols('u_i v_i')

# Construct the right-hand side of the derived equation
rhs = beta * u_i * (1 - ((1 + rho) * (1 - eta)) / (1 + phi * v_i))

# Create the full equation object
final_equation = sympy.Eq(tau_w * w_i_dot, rhs)

# Print the definitions of the simplified variables and the final equation
print("The steady-state analysis simplifies the dynamics of synaptic efficacy w_i. The new variables are defined as follows:")
print("\n- Presynaptic accumulator, v_i: Represents the input rate at synapse i.")
print(f"- Postsynaptic accumulator, u_i: Represents the shared postsynaptic calcium level, Y = sum over all synapses j of (w_j * v_j).")
print(f"- Constant rho: Represents the balance between LTD and LTP strengths, defined as rho = -alpha / beta.")
print("\nBased on these definitions, the derived equation for the change in synaptic efficacy is:")

# Use sympy's pretty print feature to display the equation
sympy.pprint(final_equation, use_unicode=True)

print(f"\n<<<tau_w * \\dot{{w}}_i = \\beta u_i (1 - ( (1+\\rho)*(1-\\eta) ) / (1 + \\phi v_i))>>>")
import sympy

# Define the symbols
tau_w, w_i, beta, u_i, rho, v_i = sympy.symbols('tau_w w_i beta u_i rho v_i')
t = sympy.symbols('t')

# Define the derivative of w_i with respect to t
dw_i_dt = sympy.Derivative(w_i, t)

# Build the equation based on the derivation
# tau_W * dW_i/dt = beta * (u_i - rho * v_i)
lhs = tau_w * dw_i_dt
rhs = beta * (u_i - rho * v_i)
equation = sympy.Eq(lhs, rhs)

# Print the final derived expression
# The problem is to derive an expression for tau_w * dw_i/dt.
# The derived expression is beta * (u_i - rho * v_i)
# The full equation is tau_w * dw_i/dt = beta * (u_i - rho * v_i)
# We will print the equation with each symbol.
print(f"The derived expression for the dynamics of synaptic efficacy is:")
print(f"{sympy.pretty(equation.lhs)} = {beta} * ({u_i} - {rho} * {v_i})")

# The final answer format as requested.
final_answer = "tau_w*d(w_i)/dt = beta*(u_i - rho*v_i)"
print(f'<<<"{final_answer}">>>')

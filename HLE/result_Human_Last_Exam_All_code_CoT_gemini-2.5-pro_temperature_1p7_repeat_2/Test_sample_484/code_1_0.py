import sympy

# Define the symbols
tau_w, w_i, t, u_i, beta, rho, v_i = sympy.symbols('tau_w w_i t u_i beta rho v_i')
d_wi_dt = sympy.Derivative(w_i, t)

# Construct the equation
lhs = tau_w * d_wi_dt
rhs = u_i * (beta + rho / (1 + v_i))
equation = sympy.Eq(lhs, rhs)

# Print the final equation
# The problem asks to print each part of the final equation.
print('The final expression is:')
print('tau_w * dw_i/dt', '=', 'u_i', '*', '(', 'beta', '+', 'rho', '/', '(', '1', '+', 'v_i', ')', ')')

# For a more readable format, uncomment the line below
# print(sympy.pretty(equation))

final_expression = "tau_w*d(w_i)/dt = u_i*(beta + rho/(1 + v_i))"
# The required output format is the expression itself in the angle brackets.
# The right hand side of the equation represents the derived expression for tau_w*dw_i/dt
print(f"\n<<<{u_i * (beta + rho / (1 + v_i))}>>>")
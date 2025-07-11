import sympy as sp

# Define the symbolic variables
q, T = sp.symbols('q T')
n = sp.Symbol('n', integer=True)

# 1. Define the general symbolic expression for the q-derivative of T^n
# The q-analogue of n, [n]_q, is (q**n - 1)/(q - 1)
q_n_symbolic = (q**n - 1) / (q - 1)
result_symbolic = q_n_symbolic * T**(n - 1)

# Print the general formula
print("The general symbolic formula for the q-derivative of T^n is:")
# Pretty print the equation
lhs = "nabla_q(T^n)"
rhs_coeff_str = sp.pretty((q**n - 1) / (q - 1), use_unicode=True)
rhs_T_str = "T^(n-1)"
print(f"{lhs} = ({rhs_coeff_str}) * {rhs_T_str}")

# As requested, output each symbolic "number" in the general equation
print("\nThe symbolic 'numbers' (components) in this equation are:")
print(f"The original power of T: {n}")
print(f"The new power of T: {n-1}")

# 2. Provide a concrete example for a specific integer, let's say n=4
print("\n" + "="*50)
print("Example for the specific integer n = 4:")
n_val = 4
q_n_val = (q**n_val - 1) / (q - 1)
# Expand the polynomial for a nice-looking output
q_n_val_expanded = sp.expand(q_n_val)
result_val_expr = q_n_val_expanded * T**(n_val - 1)

# Print the final equation for n=4
lhs_val = f"nabla_q(T^{n_val})"
rhs_val_coeff_str = sp.pretty(q_n_val_expanded, use_unicode=True)
rhs_val_T_str = f"T^({n_val - 1})"
print(f"{lhs_val} = ({rhs_val_coeff_str}) * {rhs_val_T_str}")

# Output each number in the specific equation for n=4
print("\nThe numbers in this specific equation are:")
print(f"The original power of T: n = {n_val}")
coeffs = sp.Poly(q_n_val_expanded, q).all_coeffs()
print(f"The coefficients of the q-analogue [{n_val}]_q: {[int(c) for c in coeffs]}")
print(f"The new power of T: n-1 = {n_val - 1}")

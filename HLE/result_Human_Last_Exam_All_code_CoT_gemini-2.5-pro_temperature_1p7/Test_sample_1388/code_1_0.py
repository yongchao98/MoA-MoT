import sympy

# The analysis shows that H(t) is a constant.
H_t_value = 1

# Define symbols for a nice display
u, t, u0, L2 = sympy.symbols("u, t, u_0, L^2")
norm_u = sympy.Function(f"||{u}(·,{t})||")(L2)
norm_u0 = sympy.Function(f"||{u0}||")(L2)
H = sympy.Function("H")(t)

# Formulate the inequality
inequality_general = sympy.Le(norm_u, norm_u0 * H)
inequality_specific = sympy.Le(norm_u, norm_u0 * H_t_value)

print("The explicit form of H(t) is a constant function:")
print(f"H(t) = {H_t_value}")
print("\nThe final inequality is:")
print(f"{sympy.pretty(inequality_specific, use_unicode=False)}")

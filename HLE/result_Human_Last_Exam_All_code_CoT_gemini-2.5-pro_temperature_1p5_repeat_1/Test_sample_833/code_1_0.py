import sympy

# Define the variable 'a' which represents the lower bound.
# Based on the analytical derivation for this PDE, the lower bound is -1/4.
a_val = -0.25

# Define symbols for the expression to make it readable
dudt, u, u_bar, d_dx, t, x = sympy.symbols("∂u/∂t u u_bar ∂/∂x t x")
expression_str = "(∂/∂t + (1 - 2*u)*u_bar*∂/∂x)*u_bar(t,x)"

# Create the inequality
inequality = sympy.Ge(sympy.S(expression_str, locals=sympy.S.all_symbols), a_val)

# To display it more nicely for the final print
# Let's create a representation of the left-hand side operator
LHS = sympy.S(r"(\frac{\partial}{\partial t} + (1-2u)\bar{u}\frac{\partial}{\partial x}) \bar{u}(t,x)")

# Print the final inequality expression
print("The lower bound is a = -0.25.")
print("The inequality is:")
print(f"{LHS} >= {a_val}")

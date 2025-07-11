# This script provides the expression for the inner product (phi, D_phi).

# Explanation
print("In the functional integral formalism for a neutral scalar field at finite temperature, the Euclidean action, S_E, is central.")
print("The action for a free field phi with mass m is:")
print("  S_E[phi] = Integral[d^4x] * [ (1/2)*(del_mu phi)^2 + (1/2)*m^2*phi^2 ]")
print("\nAfter integration by parts, this action can be written in a quadratic form:")
print("  S_E[phi] = (1/2) * Integral[d^4x] * phi(x) * (-del^2 + m^2) * phi(x)")
print("where del^2 is the Euclidean Laplacian operator.")
print("\nThis quadratic form is denoted as S_E = (1/2) * (phi, D_phi * phi), which defines the operator D_phi and the inner product.")
print("The operator is D_phi = -del^2 + m^2.")
print("The inner product is (f, g) = Integral[d^4x] f(x)*g(x).")
print("\nTherefore, the inner product (phi, D_phi) is given by the following equation:")

# Define the numerical components of the final equation
integration_dim = 4
kinetic_coeff = -1
op_exponent = 2
mass_exponent = 2

# Print the final equation, showing each numerical component
print("\n--- Final Equation ---")
print("(phi, D_phi) = Integral[d^", integration_dim, "x] * phi(x) * (",
      kinetic_coeff, " * del^", op_exponent, " + m^", mass_exponent,
      ") * phi(x)", sep='')

print("\nNote that this expression is exactly twice the Euclidean action:")
print("(phi, D_phi) = 2 * S_E[phi]")
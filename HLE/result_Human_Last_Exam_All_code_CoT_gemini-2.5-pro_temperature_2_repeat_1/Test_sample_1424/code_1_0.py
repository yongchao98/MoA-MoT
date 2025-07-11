# This script prints the formula for the bare Green's function G_0,
# explicitly showing its dependence on the single-particle energy epsilon_k.

print("The functional dependence of the bare Green's function G_0 on the single-particle energy epsilon_k is an inverse linear relationship.")
print("The general formula is:")
print("\n  G_0(k, omega) = 1 / (omega - (epsilon_k - mu) + i*delta)\n")
print("Where:")
print("  G_0      : The bare Green's function.")
print("  k        : The quantum number of the single-particle state (e.g., momentum).")
print("  omega    : The frequency.")
print("  epsilon_k: The single-particle energy eigenvalue for state k.")
print("  mu       : The chemical potential.")
print("  i        : The imaginary unit.")
print("  delta    : A small infinitesimal whose sign depends on the state's energy relative to the chemical potential.")
print("\nThis general form can be broken down into two specific cases:")

# Case 1: Particle state (unoccupied level)
print("\n1. For a particle state where epsilon_k > mu (an unoccupied level):")
print("   The term delta is positive (+eta, where eta is a small positive number).")
print("   The final equation is:")
print("   G_0 = 1 / (omega - epsilon_k + mu + i*eta)")

# Case 2: Hole state (occupied level)
print("\n2. For a hole state where epsilon_k < mu (an occupied level):")
print("   The term delta is negative (-eta).")
print("   The final equation is:")
print("   G_0 = 1 / (omega - epsilon_k + mu - i*eta)")

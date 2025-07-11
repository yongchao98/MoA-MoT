import numpy as np

# This script provides a simple numerical illustration of the final step of the proof.
# The proof establishes that if a sub-gradient 'phi' and a super-gradient 'psi'
# both exist at a point, they must be the same element.

# Let's represent these gradients as simple vectors for this example.
# We choose a vector for phi, an element of the sub-differential.
phi = np.array([12.5, -4.0, 7.8])

# According to the proof, if the super-differential is also non-empty, any element
# 'psi' from it must be identical to 'phi'.
psi = np.array([12.5, -4.0, 7.8])

print("--- Illustration of the Proof's Conclusion ---")
print(f"Let's assume a sub-gradient phi exists, for example: phi = {phi}")
print(f"Let's assume a super-gradient psi exists.")
print("The rigorous proof concludes that psi must be equal to phi.")
print(f"Therefore, psi must also be: psi = {psi}\n")

# The final, critical equation in the proof is that the squared norm of the
# difference vector w = psi - phi must be zero.
# ||psi - phi||^2 = 0
w = psi - phi

# In vector components, this equation is:
# (psi_1 - phi_1)^2 + (psi_2 - phi_2)^2 + ... = 0

print("The key equation derived in the proof is: ||psi - phi||^2 = 0")
print("This can be written as a sum of squares of the differences in each component.")
print("For our chosen vectors, this equation is:")

# We will build and print the equation with the specific numbers.
equation_lhs_terms = []
for i in range(len(phi)):
    term = f"({psi[i]} - {phi[i]})^2"
    equation_lhs_terms.append(term)
final_equation_str = " + ".join(equation_lhs_terms) + " = 0"
print(final_equation_str)

print("\nLet's evaluate the left side of the equation:")
evaluated_terms = []
for i in range(len(w)):
    evaluated_terms.append(f"({w[i]})^2")
calculation_str = " + ".join(evaluated_terms)
result = np.dot(w, w)
print(f"  {calculation_str} = {result}")

print("\nBecause a sum of non-negative square terms is zero, each individual term must be zero.")
for i in range(len(phi)):
    print(f"  ({psi[i]} - {phi[i]})^2 = 0  =>  {psi[i]} = {phi[i]}")

print("\nThis demonstrates that the vector psi must be component-wise identical to the vector phi.")
print("Thus, if both the sub-differential and super-differential are non-empty, they must be the same singleton set, which implies the functional is differentiable.")
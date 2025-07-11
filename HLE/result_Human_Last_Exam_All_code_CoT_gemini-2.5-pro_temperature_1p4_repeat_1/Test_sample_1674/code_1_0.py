import numpy as np

# Set a seed for reproducibility of the 'random' medium
np.random.seed(0)

# 1. Define the optical elements and the input beam

# Input beam with horizontal polarization
E_in = np.array([[1], [0]], dtype=complex)

# A random medium described by a 2x2 complex transmission matrix 'T'
# This represents a medium that arbitrarily alters polarization.
T = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)

# The inverse of the random medium's matrix 'T'.
# The initial theory assumes we can use this to reverse the effect of T.
try:
    T_inv = np.linalg.inv(T)
except np.linalg.LinAlgError:
    print("The random matrix T is not invertible. Rerunning.")
    exit()

# A birefringent medium. Let's use a quarter-wave plate with its fast axis vertical.
# This introduces a pi/2 phase shift between horizontal and vertical components.
J_b = np.array([[1, 0], [0, 1j]], dtype=complex)

# 2. Simulate the process
# The user's question is whether the theory holds if we add a birefringent layer.
# Let's model the new system where the light first passes through the random medium 'T'
# and then through the birefringent plate 'J_b'.
# The output of this combined system is: E_out = J_b @ T @ E_in

# Now, we test the "theory" by applying the inverse of the original random medium, T_inv.
# According to the theory, this should recover the input needed to produce the output.
# Let's see what happens.
E_final = T_inv @ (J_b @ T @ E_in)

# 3. Print the results and explain
print("This script simulates an optical system to test if a theory of inversion holds.")
print("The theory suggests that the effect of a random medium (T) can be undone by its inverse (T_inv).")
print("We test if this holds when a birefringent plate (J_b) is added to the system.")
print("-" * 70)

print("The equation we are testing is: T_inv * (J_b * T * E_in) = E_in ?")
print("-" * 70)
print("The components of the equation are:")

print("\nOriginal Input Polarization Vector (E_in):")
print(E_in)

print("\nRandom Medium Matrix (T):")
print(T)

print("\nBirefringent Plate Matrix (J_b, a quarter-wave plate):")
print(J_b)

print("\nInverse of Random Medium Matrix (T_inv):")
print(T_inv)
print("-" * 70)

print("Let's calculate the final vector E_final = T_inv @ J_b @ T @ E_in")
print(f"E_final = \n{T_inv} \n\n@ \n\n{J_b} \n\n@ \n\n{T} \n\n@ \n\n{E_in}\n\n= \n{E_final}")
print("-" * 70)

# 4. Final Conclusion
print("Conclusion:")
print("\nOriginal Input E_in:\n", E_in)
print("\nFinal Calculated Vector E_final:\n", E_final)

# Compare the final result with the initial input
are_close = np.allclose(E_in, E_final)

print(f"\nIs the final vector equal to the original input? {are_close}")

if not are_close:
    print("\nAs shown, the final vector is NOT the same as the original input.")
    print("This demonstrates that the theory does not hold.")
    print("The reason is that matrix multiplication is not commutative. The term (T_inv @ J_b @ T) does not simplify to the identity matrix.")
    print("Therefore, simply applying the inverse of the random medium is not sufficient to reverse the effects of the combined system.")
else:
    # This case is extremely unlikely for random matrices
    print("\nIn this highly unlikely case, the vectors are the same.")

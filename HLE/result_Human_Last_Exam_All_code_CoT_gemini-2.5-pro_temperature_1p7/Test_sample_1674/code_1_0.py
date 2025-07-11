import numpy as np

def print_vector(name, vec):
    """Helper function to print a complex vector."""
    print(f"{name}:\n[ {vec[0]:.2f} ]\n[ {vec[1]:.2f} ]\n")

def print_matrix(name, mat):
    """Helper function to print a complex matrix."""
    print(f"{name}:\n[ {mat[0,0]:.2f}  {mat[0,1]:.2f} ]\n[ {mat[1,0]:.2f}  {mat[1,1]:.2f} ]\n")

# --- Setup ---
# Fix the random seed for reproducibility
np.random.seed(42)

# 1. Define the input beam: A Laguerre-Gaussian beam with horizontal polarization.
# We represent its polarization as a Jones vector.
E_in = np.array([1 + 0j, 0 + 0j])
print("--- Initial Setup ---")
print_vector("Input Beam Polarization (E_in)", E_in)

# 2. Define the random medium's transmission matrix (T_rand).
# We create a random 2x2 complex matrix.
T_rand = np.random.rand(2, 2) + 1j * np.random.rand(2, 2)
print_matrix("Random Medium TM (T_rand)", T_rand)

# 3. Define the birefringent medium's matrix (B).
# This is a quarter-wave plate with its fast axis vertical.
# It induces a pi/2 (or i) phase shift on the vertical component.
B = np.array([[1 + 0j, 0 + 0j],
              [0 + 0j, 0 + 1j]])
print_matrix("Birefringent Medium (B)", B)


# --- Scenario 1: No Birefringence ---
print("\n--- Scenario 1: System with Random Medium ONLY ---")
# Propagate through the random medium
E_out_1 = T_rand @ E_in
# Find the inverse of the random medium's matrix
T_rand_inv = np.linalg.inv(T_rand)
# Propagate the output back through the inverse matrix
E_recovered = T_rand_inv @ E_out_1

print("Verifying recovery: T_rand_inv * (T_rand * E_in) = E_in")
print("Result:")
print_vector("Recovered Beam (E_recovered)", E_recovered)
print("The recovered beam matches the input beam, as expected.")


# --- Scenario 2: With Birefringence ---
print("\n--- Scenario 2: System with Birefringent + Random Medium ---")
# The theory states that to recover the input, one must use the inverse of the
# random medium's TM. Let's test if this theory holds when we add a birefringent layer.

# The full system matrix is now T_rand * B
T_total = T_rand @ B

# The output of the full system is T_total * E_in
E_out_2 = T_total @ E_in

# Now, we apply the inverse of ONLY the random medium part, as per the original theory.
E_failed_recovery = T_rand_inv @ E_out_2

print("Testing the theory: T_rand_inv * (T_total * E_in) should equal E_in?")
print("This is equivalent to: T_rand_inv * (T_rand * B * E_in) = E_in")
print("Let's calculate the final resulting vector of this operation.")

print("\nFinal Equation and Calculation:")

print_matrix("T_rand_inv", T_rand_inv)
print("     * ")
print_matrix("T_total = T_rand * B", T_total)
print("     * ")
print_vector("E_in", E_in)
print("     = ")
print_vector("Final Result", E_failed_recovery)

print("\n--- Conclusion ---")
print("As you can see, the final vector is not the original input [1, 0].")
print("Therefore, the theory that you can reverse the effects of the combined system by simply using the inverse of the random medium's transmission matrix does not hold.")
print("The birefringent element alters the polarization state before it enters the random medium, and this effect is not accounted for in the simple reversal step.")

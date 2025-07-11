import numpy as np

def generate_random_unitary():
    """Generates a random 2x2 unitary matrix T."""
    # Create a random complex matrix Z
    Z = np.random.randn(2, 2) + 1j * np.random.randn(2, 2)
    # Compute the QR decomposition of Z
    Q, R = np.linalg.qr(Z)
    return Q

# --- Step 1: Define the optical components and input beam ---

# Let's assume a horizontally polarized input beam.
# This is our Jones vector for the input polarization.
p_in = np.array([1, 0])

# T: A random Jones matrix for the random medium.
# We generate a random unitary matrix to represent a lossless random medium.
T = generate_random_unitary()
T_inv = np.linalg.inv(T) # The inverse of the random medium

# B: A Jones matrix for a birefringent medium.
# Let's model a standard quarter-wave plate with its fast axis horizontal.
# It shifts the phase of the vertical component by 90 degrees (pi/2 radians).
B = np.array([[1, 0],
              [0, 1j]])

print("--- System Setup ---")
print(f"Input Polarization (p_in):\n{p_in}\n")
print(f"Random Medium (T):\n{np.round(T, 3)}\n")
print(f"Birefringent Plate (B):\n{B}\n")

# --- Step 2: Simulate the forward path to get output1 ---
# output1 = S' * p_in = B * T * p_in
output1 = B @ T @ p_in

print("--- Forward Propagation ---")
print("output1 = B * T * p_in")
print(f"output1:\n{np.round(output1, 3)}\n")


# --- Step 3: Simulate the described "reversal" to get output2 ---
# According to the problem, this reversal only inverts T.
# output2 = T_inverse * output1
output2 = T_inv @ output1

print("--- Assumed 'Reversal' to find necessary input ---")
print("output2 = T_inverse * output1")
print(f"output2 (the proposed necessary input):\n{np.round(output2, 3)}\n")

# --- Step 4: Test the theory ---
# The theory states that `output2` should be the input that gives `output1`.
# We test this by calculating `test_output = S' * output2 = B * T * output2`.
# If the theory holds, test_output must equal output1.

test_output = B @ T @ output2

# As shown in the plan, this is equivalent to checking if B^2 * T * p_in == B * T * p_in
LHS = B @ B @ T @ p_in
RHS = B @ T @ p_in

print("--- Testing the Theory ---")
print("The theory holds if (B * T * output2) equals output1.")
print("This is mathematically equivalent to checking if the equation 'B*B*T*p_in = B*T*p_in' is true.")
print("\nCalculating the Left Hand Side (LHS = B*B*T*p_in):")
print(f"LHS = \n{np.round(LHS, 3)}")

print("\nCalculating the Right Hand Side (RHS = B*T*p_in), which is just output1:")
print(f"RHS = \n{np.round(RHS, 3)}")

# --- Step 5: Conclusion ---
is_equal = np.allclose(LHS, RHS)
print("\n--- Conclusion ---")
print(f"Are the two sides of the equation equal? {is_equal}")

if not is_equal:
    print("The results are different. This demonstrates that when the birefringent medium is added,")
    print("the described reversal procedure fails. The theory does not hold because the reversal")
    print("procedure does not account for the effect of the birefringent plate.")
else:
    print("The results are the same. This would only happen in trivial cases (e.g., B is the identity matrix).")

<<<Yes, the theory can fail.>>>
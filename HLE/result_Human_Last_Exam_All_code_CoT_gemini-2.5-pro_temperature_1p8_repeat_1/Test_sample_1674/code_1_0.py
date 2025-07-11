import numpy as np

def print_vector(name, vec):
    """Helper function to print complex vectors for clarity."""
    # Rounding to 3 decimal places for readability
    rounded_vec = np.round(vec, 3)
    print(f"{name}:")
    print(rounded_vec)
    print("-" * 25)

# --- General Setup ---
# Define the "space" of our beam with a set number of modes (like pixels)
N_MODES = 4
# Use a fixed seed for random numbers to get reproducible results
np.random.seed(0)

# Define the initial input beam.
# It's purely horizontally polarized and exists only in the first spatial mode.
# In our vector format, the first N_MODES elements are the X-polarization component
# and the next N_MODES elements are the Y-polarization component.
input_beam_x = np.zeros(N_MODES, dtype=complex)
input_beam_x[0] = 1.0
input_beam_y = np.zeros(N_MODES, dtype=complex)
full_input_beam = np.concatenate((input_beam_x, input_beam_y))

print("--- Initial Setup ---")
print_vector("Initial Input (X-pol)", input_beam_x)


# --- Case 1: Scalar System (The Theory Holds) ---
print("\n--- Case 1: Scalar System (No Birefringence) ---")
print("Propagating through a random medium. Polarization is unaffected.")

# Create a random scalar transmission matrix for the scattering medium.
T_scalar = np.random.randn(N_MODES, N_MODES) + 1j * np.random.randn(N_MODES, N_MODES)

# The final equation for this step is: output = T_scalar * input
output_beam_x = T_scalar @ input_beam_x
print_vector("Scrambled Output (X-pol)", output_beam_x)

# Now, we apply the theory of inversion.
# Equation: retrieved_input = T_scalar_inverse * output
T_scalar_inv = np.linalg.inv(T_scalar)
retrieved_input_x = T_scalar_inv @ output_beam_x

print("Applying the inverse of the scalar TM to the output...")
print_vector("Retrieved Input (X-pol)", retrieved_input_x)
print("SUCCESS: The theory holds. The retrieved input matches the original.")


# --- Case 2: Vectorial System (The Theory Fails) ---
print("\n\n--- Case 2: System With Birefringence ---")
print("Adding a birefringent plate, then attempting the same SCALAR inversion.")

# Model a birefringent element (quarter-wave plate at 45 deg) that mixes polarizations.
# Jones matrix for the plate:
J_birefringent = 0.5 * np.array([[1, 1j], [-1j, 1]])

# Create the full vectorial matrices for the system components.
# The random medium doesn't mix polarization, so its matrix is block-diagonal.
T_s_vec = np.block([[T_scalar, np.zeros((N_MODES, N_MODES))],
                    [np.zeros((N_MODES, N_MODES)), T_scalar]])

# The birefringent plate acts on polarization for all spatial modes.
T_b_vec = np.kron(J_birefringent, np.eye(N_MODES))

# The full system is the random medium followed by the plate.
T_total_vec = T_b_vec @ T_s_vec

# The final equation for propagation is: full_output = T_total * full_input
full_output_beam = T_total_vec @ full_input_beam
output_beam_x_vec = full_output_beam[:N_MODES]
output_beam_y_vec = full_output_beam[N_MODES:]

print("Propagating through Random Medium + Birefringent Plate...")
print_vector("Final Output (X-pol part)", output_beam_x_vec)
print_vector("Final Output (Y-pol part)", output_beam_y_vec)
print("Note: The Y-pol component is no longer zero due to polarization scrambling.")

# THE CRITICAL TEST: Attempt to invert using the original, simple scalar theory.
# This simulates an experimenter who is unaware of the birefringent plate and
# incorrectly assumes the system is scalar.
# Equation: retrieved_input = T_scalar_inverse * output_x_part
print("\nAttempting inversion using the INCOMPLETE scalar model...")
failed_retrieval_x = T_scalar_inv @ output_beam_x_vec
print("The final equation we are testing is:")
print("retrieved_x = T_scalar_inv @ final_output_x")
print_vector("Result of Equation", failed_retrieval_x)
print("FAILURE: The theory does not hold. The retrieved input does NOT match the original.")
print("\nThis is because the scalar model (T_scalar_inv) cannot account for the\n"
      "polarization scrambling introduced by the birefringent element.")

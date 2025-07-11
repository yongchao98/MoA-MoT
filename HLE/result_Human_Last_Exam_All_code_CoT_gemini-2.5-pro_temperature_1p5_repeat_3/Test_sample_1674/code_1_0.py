import numpy as np

def pretty_print_vector(v, name="Vector"):
    """Helper function to print a complex vector nicely."""
    print(f"{name}:\n[ {v[0,0]:.2f} {v[0,0].imag:+.2f}j ]\n[ {v[1,0]:.2f} {v[1,0].imag:+.2f}j ]\n")

def pretty_print_matrix(m, name="Matrix"):
    """Helper function to print a complex matrix nicely."""
    print(f"{name}:")
    for row in m:
        print(f"[ {row[0]:.2f} {row[0].imag:+.2f}j   {row[1]:.2f} {row[1].imag:+.2f}j ]")
    print()

# 1. Define the initial state: a horizontally polarized Laguerre-Gaussian beam.
# We only model the polarization state using a Jones vector.
# [1] represents horizontal polarization, [0] represents vertical.
V_in = np.array([[1], [0]], dtype=complex)
print("--------------------------------------------------")
print("Step 1: Define the input beam's polarization.")
print("--------------------------------------------------")
pretty_print_vector(V_in, "Input Polarization Vector (V_in)")

# 2. Define the optical elements' transmission matrices.
# The random medium T is modeled as a lossless (unitary) scrambler.
# A rotation matrix is a simple example of a unitary matrix.
theta_T = np.deg2rad(15)
T = np.array([
    [np.cos(theta_T), -np.sin(theta_T)],
    [np.sin(theta_T), np.cos(theta_T)]
])

# The birefringent medium B has a phase shift (birefringence) and
# polarization-dependent loss (dichroism), making it non-unitary.
# It shifts the phase of the vertical component by 90 degrees (pi/2)
# and also absorbs 50% of its amplitude.
B = np.array([
    [1, 0],
    [0, 0.5 * np.exp(1j * np.pi / 2)] # 0.5 for loss, exp() for phase shift
])

print("--------------------------------------------------")
print("Step 2: Define the system's components.")
print("--------------------------------------------------")
pretty_print_matrix(T, "Random Medium Matrix (T)")
pretty_print_matrix(B, "Birefringent/Dichroic Plate Matrix (B)")

# 3. Calculate the full system matrix S and the output state V_out.
# The beam passes through T, then B.
S = B @ T
V_out = S @ V_in

print("--------------------------------------------------")
print("Step 3: Calculate the output of the forward pass.")
print("Equation: V_out = B * T * V_in")
print("--------------------------------------------------")
pretty_print_matrix(S, "Total System Matrix (S = B @ T)")
pretty_print_vector(V_out, "Output Polarization Vector (V_out)")

# 4. Simulate the reversal using optical phase conjugation (OPC).
# For a reciprocal medium, the reverse propagation path is described by the
# transpose of the system matrix, S_T.
# The OPC process creates a phase-conjugated copy of the output vector.
S_T = S.T
V_out_conj = np.conjugate(V_out)

print("--------------------------------------------------")
print("Step 4: Prepare for the reversal (backward pass).")
print("--------------------------------------------------")
pretty_print_matrix(S_T, "Reverse System Matrix (S_T = Transpose of S)")
pretty_print_vector(V_out_conj, "Phase-Conjugated Output Vector (V_out*)")

# 5. Calculate the final state after the reversal.
# Equation for reversal: V_final = S_T * V_out_conj
V_final = S_T @ V_out_conj

# The "theory" holds if the final state is the conjugate of the *original* input.
V_in_conj = np.conjugate(V_in)

print("--------------------------------------------------")
print("Step 5: Calculate the final vector after reversal.")
print("Equation: V_final = S_T * V_out*")
print("--------------------------------------------------")
pretty_print_vector(V_final, "Final Reversed Vector (V_final)")
pretty_print_vector(V_in_conj, "Desired Result (V_in*)")

# 6. Conclusion: Compare the result to the desired outcome.
# We use np.allclose to check if the two vectors are equal within a small tolerance.
is_reversal_successful = np.allclose(V_final, V_in_conj)

print("--------------------------------------------------")
print("Step 6: Conclusion.")
print("--------------------------------------------------")
print(f"Does the final vector match the desired (conjugated input) vector? {is_reversal_successful}")
if not is_reversal_successful:
    print("The reversal failed. The theory does not hold because the birefringent element introduced loss (dichroism), making the system non-unitary.")
else:
    print("The reversal was successful.")

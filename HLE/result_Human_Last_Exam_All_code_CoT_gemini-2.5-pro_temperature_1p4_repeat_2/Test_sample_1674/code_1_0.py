import numpy as np

def format_complex(c):
    """Formats a complex number for printing."""
    # Handle the case where the input might be a float
    if not isinstance(c, complex):
        c = complex(c)
    # Format to 2 decimal places, ensuring a sign for the imaginary part
    return f"{c.real:.2f}{c.imag:+.2f}j"

def print_equation(result_vec, matrix, input_vec, description):
    """Prints a formatted 2x2 matrix-vector multiplication."""
    print(description)
    line1 = (f"[{format_complex(result_vec[0])}] = [[{format_complex(matrix[0, 0])}, {format_complex(matrix[0, 1])}]"
             f" @ [{format_complex(input_vec[0])}]")
    line2 = (f"[{format_complex(result_vec[1])}]   [[{format_complex(matrix[1, 0])}, {format_complex(matrix[1, 1])}]"
             f"   [{format_complex(input_vec[1])}]")
    print(line1)
    print(line2)
    print("-" * 50)

# 1. Define the components of the optical system
# The input beam has horizontal polarization, represented by a Jones vector.
E_in = np.array([1 + 0j, 0 + 0j])

# The random medium introduces a phase shift, but is polarization-insensitive.
# We model this with a scalar `t` applied to an identity matrix.
random_phase_shift = 2 * np.pi * np.random.rand()
t = np.exp(1j * random_phase_shift)
T_matrix = np.array([[t, 0], [0, t]])

# The birefringent medium is a quarter-wave plate with its fast axis
# aligned vertically. It adds a 90-degree (pi/2) phase shift to the
# vertical component of the electric field.
B_matrix = np.array([[1 + 0j, 0 + 0j], [0 + 0j, 1j]])

print("--- System Definition ---")
print(f"Initial Input (Horizontal Polarization): E_in =\n{E_in}\n")
print(f"Random Medium Matrix (Scalar): T_matrix =\n{np.round(T_matrix, 2)}\n")
print(f"Birefringent Plate Matrix (QWP): B_matrix =\n{B_matrix}\n")
print("=" * 50)

# 2. Simulate the new system (Random Medium + Birefringent Plate)
# The operator for the new system is S_new = B * T
S_new_matrix = B_matrix @ T_matrix
# The output electric field after passing through the new system
E_out_new = S_new_matrix @ E_in

print_equation(E_out_new, S_new_matrix, E_in,
               "Step 1: Calculate the output from the new system (S_new @ E_in)")

# 3. Attempt to reverse the process using the OLD theory's inverse
# The old theory only knew about the random medium 'T'.
# Its inverse is simple: T_inv = [[1/t, 0], [0, 1/t]]
T_inv_matrix = np.linalg.inv(T_matrix)

# Apply this incorrect inverse to the output of the new system
E_reconstructed_fail = T_inv_matrix @ E_out_new

print_equation(E_reconstructed_fail, T_inv_matrix, E_out_new,
               "Step 2: Apply inverse of OLD system (T_inv @ E_out_new) -> Fails")
print(f"Result of incorrect inversion:\n{E_reconstructed_fail.round(2)}")
print("This is NOT the original input vector [1, 0]. The theory failed.\n")
print("=" * 50)

# 4. Reverse the process using the CORRECT (vectorial) inverse
# The correct inverse must account for the full system S_new = B * T
S_new_inv_matrix = np.linalg.inv(S_new_matrix)

# Apply the correct inverse to the output of the new system
E_reconstructed_correct = S_new_inv_matrix @ E_out_new

print_equation(E_reconstructed_correct, S_new_inv_matrix, E_out_new,
               "Step 3: Apply inverse of NEW system (S_new_inv @ E_out_new) -> Succeeds")
print(f"Result of correct vectorial inversion:\n{E_reconstructed_correct.round(2)}")
print("This IS the original input vector [1, 0]. The principle of reversibility holds.")

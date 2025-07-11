import numpy as np

def print_results(device_name, E_in, E_final):
    """Helper function to print the results in a readable format."""
    print(f"--- Testing: {device_name} ---")
    print(f"Initial Polarization: [{E_in[0,0]:.2f}, {E_in[1,0]:.2f}]")
    print(f"Final Polarization (after round trip): [{E_final[0,0]:.2f}, {E_final[1,0]:.2f}]")
    # np.allclose checks if two arrays are element-wise equal within a tolerance.
    if np.allclose(E_in, E_final):
        print("Result: Initial polarization was recovered. The theory holds.\n")
    else:
        print("Result: Initial polarization was NOT recovered. The theory fails.\n")

# --- Setup ---
# Define an initial horizontal polarization state E = [Ex, Ey]
# Using column vectors for matrix multiplication: [[1], [0]]
E_initial = np.array([[1], [0]])

# Angle for the devices (45 degrees = pi/4 radians)
theta = np.pi / 4

# --- Case 1: Reciprocal Device (Half-Wave Plate at 45 degrees) ---
# A HWP at 45 degrees swaps horizontal and vertical polarizations.
# The forward matrix T_fwd is: [[cos(2*theta), sin(2*theta)], [sin(2*theta), -cos(2*theta)]]
c2t = np.cos(2 * theta)
s2t = np.sin(2 * theta)
T_hwp_fwd = np.array([[c2t, s2t], [s2t, -c2t]])

# For a reciprocal device, the reverse path matrix is the transpose of the forward matrix.
T_hwp_rev = T_hwp_fwd.T

# Simulate the round trip: pass through forward, then backward.
E_out_hwp = T_hwp_fwd @ E_initial
E_final_hwp = T_hwp_rev @ E_out_hwp

# Print results for the HWP
print_results("Reciprocal Half-Wave Plate", E_initial, E_final_hwp)


# --- Case 2: Non-Reciprocal Device (Faraday Rotator with 45 degree rotation) ---
# The forward matrix for a rotation of 'theta' is T_fwd:
ct = np.cos(theta)
st = np.sin(theta)
T_fr_fwd = np.array([[ct, -st], [st, ct]])

# For a non-reciprocal Faraday rotator, the reverse path matrix T_rev is the SAME
# as the forward matrix. The rotation is in the same direction regardless of travel direction.
T_fr_rev = T_fr_fwd

# Simulate the round trip: pass through forward, then backward.
E_out_fr = T_fr_fwd @ E_initial
E_final_fr = T_fr_rev @ E_out_fr

# Print results for the Faraday Rotator
print_results("Non-Reciprocal Faraday Rotator", E_initial, E_final_fr)
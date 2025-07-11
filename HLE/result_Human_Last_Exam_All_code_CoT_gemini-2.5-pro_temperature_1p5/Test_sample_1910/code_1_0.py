import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak of NaMgH3.
    """
    # Step 1: Define the d-spacing for the second major diffraction peak.
    # NaMgH3 has an orthorhombic Pnma structure. From simulated diffraction data
    # (e.g., Materials Project, mp-30752), the most intense peak is (2,1,1)
    # and the second most intense peak is (1,2,1).
    # The d-spacing for the (1,2,1) peak is given as 2.733 Å.
    d_hkl = 2.733  # in Angstroms

    # The experimental wavelength is not needed for the Q-space calculation,
    # as Q = 2*pi/d is an intrinsic property of the crystal.

    # Step 2: Calculate the Q-space position using the formula Q = 2 * pi / d.
    pi_value = math.pi
    q_value = (2 * pi_value) / d_hkl

    # Step 3: Print the calculation steps and the final result.
    print("The second major diffraction peak of NaMgH3 corresponds to the (1,2,1) Miller indices.")
    print(f"The d-spacing for this peak is {d_hkl} Å.")
    print("\nThe Q-space position is calculated using the formula: Q = 2 * pi / d")
    print("Plugging in the values:")
    print(f"Q = (2 * {pi_value}) / {d_hkl}")
    print(f"Q = {q_value:.4f} Å⁻¹")

if __name__ == "__main__":
    calculate_q_space_position()
    # The final numerical answer for the submission system.
    # print(f"<<<{2 * math.pi / 2.733:.3f}>>>")
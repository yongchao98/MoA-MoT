import numpy as np

def analyze_vector_beam_generation():
    """
    Demonstrates the impossibility of generating an arbitrary vector beam from
    a linearly polarized input through a fixed linear optical system.
    """
    # Let's model the entire complex optical system's effect at a single point
    # or for a single spatial mode using a 2x2 Jones Matrix.
    # A random medium that depolarizes light will have non-zero off-diagonal elements,
    # enabling it to convert horizontal polarization to vertical and vice-versa.
    # E_out = System_Matrix * E_in
    # |Ex_out| = |S_xx S_xy| * |Ex_in|
    # |Ey_out|   |S_yx S_yy|   |Ey_in|
    #
    # We define a fixed, arbitrary system matrix.
    system_jones_matrix = np.array([
        [0.5 + 0.1j, 0.3 - 0.2j],
        [0.4 + 0.4j, 0.2 - 0.3j]
    ])

    # The input beam has a fixed linear polarization (e.g., horizontal).
    # We can control its overall complex amplitude 'a' (phase and amplitude).
    # This means the input Jones vector is always in the form: E_in = [a, 0]^T
    #
    # The output is calculated as: E_out = system_jones_matrix @ [a, 0]^T
    # This gives:
    # Ex_out = S_xx * a
    # Ey_out = S_yx * a
    #
    # The crucial insight is the ratio of the output components:
    # Ex_out / Ey_out = (S_xx * a) / (S_yx * a) = S_xx / S_yx
    # This ratio is a constant value determined entirely by the system matrix.
    # It does not depend on our choice of the input amplitude 'a'.
    # Therefore, we can only produce output beams where the polarization state
    # (defined by the ratio of Ex to Ey) is equal to this fixed value.

    s_xx = system_jones_matrix[0, 0]
    s_yx = system_jones_matrix[1, 0]

    # This check ensures the system can create some y-component from an x-component input.
    if np.isclose(s_yx, 0):
        print("The system cannot convert x-polarization to y-polarization. Vector beam generation is impossible.")
        return

    system_ratio = s_xx / s_yx

    print(f"The optical system is represented by the Jones Matrix:\n{system_jones_matrix}")
    print("\nFor a linearly polarized input, the system can only produce an output with a fixed polarization state.")
    print(f"The equation for the ratio of output components is: Ex_out / Ey_out = S_xx / S_yx")
    print(f"For this system, the fixed complex ratio Ex_out / Ey_out is: ({s_xx}) / ({s_yx}) = {system_ratio:.4f}\n")

    print("--- Can we create an ARBITRARY vector beam? ---")

    # Test Case 1: Try to generate right-hand circularly polarized light.
    # Desired E_out = [1, 1j]^T. The required ratio is Ex/Ey = 1 / 1j = -1j.
    target1_ratio = 1 / (1j)
    print("\nTarget 1: Right-hand circular polarization (Ex/Ey = -1.0000j)")
    print(f"The system's fixed ratio is {system_ratio:.4f}. The target ratio is {target1_ratio:.4f}.")
    print("Result: This is NOT achievable because the target ratio does not match the system's fixed ratio.")

    # Test Case 2: Try to generate linearly polarized light at 45 degrees.
    # Desired E_out = [1, 1]^T. The required ratio is Ex/Ey = 1 / 1 = 1.
    target2_ratio = 1.0
    print("\nTarget 2: 45-degree linear polarization (Ex/Ey = 1.0000)")
    print(f"The system's fixed ratio is {system_ratio:.4f}. The target ratio is {target2_ratio:.4f}.")
    print("Result: This is NOT achievable because the target ratio does not match the system's fixed ratio.")

    print("\nConclusion: Since we can only produce one specific output polarization state (one ratio of components), we cannot produce an 'arbitrary' vector beam, which would require the ability to create any polarization state (any ratio).")


if __name__ == '__main__':
    analyze_vector_beam_generation()
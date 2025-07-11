import math

def calculate_q_space_peak():
    """
    Calculates the Q-space position for the second major diffraction peak of NaMgH3.
    """
    # 1. Define crystal structure parameters
    # Pseudo-cubic lattice parameter for NaMgH3 in Angstroms
    a = 3.915
    # Miller indices for the second major peak in a simple perovskite structure
    h, k, l = 2, 0, 0

    # 2. Calculate the d-spacing for the (200) peak
    # d_hkl = a / sqrt(h^2 + k^2 + l^2)
    d_spacing = a / math.sqrt(h**2 + k**2 + l**2)

    # 3. Calculate the Q-space position
    # Q = 2 * pi / d
    q_value = 2 * math.pi / d_spacing

    # 4. Print the results
    print(f"The second major diffraction peak for NaMgH3 corresponds to the ({h},{k},{l}) reflection.")
    print(f"Using a pseudo-cubic lattice parameter a = {a} Å.")
    print(f"The calculated Q-space position is: {q_value:.4f} 1/Å.")
    
    # As requested, print the final equation with all numbers
    print("\nCalculation breakdown:")
    print(f"d_spacing = {a} / sqrt({h}^2 + {k}^2 + {l}^2) = {d_spacing:.4f} Å")
    print(f"Q = (2 * {math.pi:.6f}) / {d_spacing:.4f} Å = {q_value:.4f} 1/Å")


if __name__ == "__main__":
    calculate_q_space_peak()
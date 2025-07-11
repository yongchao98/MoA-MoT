import math

def calculate_q_space_position():
    """
    Calculates the Q-space position for the second major diffraction peak
    of orthorhombic NaMgH3.
    """
    # Step 1: Define crystal structure parameters for NaMgH3 (Pnma)
    # Lattice parameters in Angstroms (Å)
    a = 5.46
    b = 7.77
    c = 5.45

    # Step 2: Identify Miller indices for the second major peak
    # This corresponds to the pseudo-cubic {200} reflection, which is
    # represented by the (202) peak in the orthorhombic cell.
    h, k, l = 2, 0, 2

    # Step 3: Calculate the interplanar spacing (d)
    # Formula for an orthorhombic crystal: 1/d^2 = (h/a)^2 + (k/b)^2 + (l/c)^2
    try:
        d_squared_inv = (h / a)**2 + (k / b)**2 + (l / c)**2
        d_spacing = 1 / math.sqrt(d_squared_inv)
    except (ValueError, ZeroDivisionError) as e:
        print(f"Error calculating d-spacing: {e}")
        return

    # Step 4: Calculate the Q-space position
    # Formula: Q = 2 * pi / d
    q_position = 2 * math.pi / d_spacing

    # Print the results, including the full equation
    print(f"The second major diffraction peak for NaMgH3 corresponds to the ({h}{k}{l}) reflection.")
    print(f"Lattice parameters used: a={a} Å, b={b} Å, c={c} Å.")
    print("\nThe Q-space position (Q) is calculated as follows:")
    print(f"1. Interplanar spacing (d) for ({h}{k}{l}): d = 1 / sqrt(({h}/{a})² + ({k}/{b})² + ({l}/{c})²) = {d_spacing:.4f} Å")
    print(f"2. Q-space position: Q = 2 * \u03C0 / d = 2 * {math.pi:.5f} / {d_spacing:.4f} = {q_position:.3f} 1/Å")
    print("\nNote: The measurement wavelength of 0.2952 Å is not needed for this calculation as Q-space is an intrinsic property of the crystal lattice.")

if __name__ == '__main__':
    calculate_q_space_position()

import math

def calculate_q_space_position():
    """
    Calculates the Q-space position of the second major diffraction peak
    of NaMgH3 perovskite.
    """
    # Step 1: Define known lattice parameters for orthorhombic NaMgH3
    a_ortho = 5.518  # Angstrom
    b_ortho = 5.632  # Angstrom
    c_ortho = 7.892  # Angstrom
    Z = 4  # Formula units per orthorhombic cell

    # Step 2: Calculate the pseudocubic lattice parameter
    # Volume of the orthorhombic cell
    V_ortho = a_ortho * b_ortho * c_ortho
    # Effective volume of one cubic formula unit
    V_cubic = V_ortho / Z
    # Pseudocubic lattice parameter
    a_pc = V_cubic**(1/3)

    # Step 3: Miller indices for the second major peak in the pseudocubic cell
    # This corresponds to the {200} family of reflections.
    h, k, l = 2, 0, 0

    # Step 4: Calculate the d-spacing for the (200) peak
    d_200 = a_pc / math.sqrt(h**2 + k**2 + l**2)

    # Step 5: Calculate the Q-space position
    Q = 2 * math.pi / d_200

    print("Step 1: The orthorhombic lattice parameters are a = {} Å, b = {} Å, c = {} Å.".format(a_ortho, b_ortho, c_ortho))
    print("Step 2: The effective pseudocubic lattice parameter 'a_pc' is calculated to be {:.4f} Å.".format(a_pc))
    print("Step 3: The second major peak corresponds to the pseudocubic (200) reflection.")
    print("Step 4: The d-spacing (d_200) for this peak is calculated to be {:.4f} Å.".format(d_200))
    print("\nFinal Calculation:")
    print("The Q-space position is given by Q = 2 * pi / d_200")
    print("Q = 2 * {:.4f} / {:.4f} Å".format(math.pi, d_200))
    print("The position of the second major diffraction peak is Q = {:.4f} 1/Å.".format(Q))


calculate_q_space_position()
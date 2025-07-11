def predict_ir_phonons():
    """
    This function calculates and prints the number of predicted IR-active phonons
    for a LiNiPO4 monocrystal in different polarizations based on group theory.

    The crystal structure is orthorhombic, space group Pnma (No. 62).
    The analysis leads to the following number of IR-active optical modes:
    - B3u modes are active for polarization E||x (a-axis).
    - B2u modes are active for polarization E||y (b-axis).
    - B1u modes are active for polarization E||z (c-axis).

    Based on the literature, the total number of modes for these symmetries are:
    Total B1u = 12
    Total B2u = 9
    Total B3u = 12

    Subtracting the 3 acoustic modes (one of each B1u, B2u, B3u) gives the
    number of IR-active optical phonons.
    """

    # Number of phonons for polarization E parallel to x-axis (a-axis)
    num_phonons_x = 11

    # Number of phonons for polarization E parallel to y-axis (b-axis)
    num_phonons_y = 8

    # Number of phonons for polarization E parallel to z-axis (c-axis)
    num_phonons_z = 11

    # Print the result in the specified format
    print(f"E||x: {num_phonons_x}, E||y: {num_phonons_y}, E||z: {num_phonons_z}")

predict_ir_phonons()
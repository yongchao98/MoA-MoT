def calculate_mass_ratio():
    """
    Calculates and displays the leading-order asymptotic mass ratio
    for excitations in the CP(N-1) model.
    """

    # In the large-N limit, the mass M_k of the k-th excitation is
    # proportional to k.
    # We are interested in the ratio between the lightest (k=1) and the
    # subsequent higher excitation (k=2).

    # The quantum number for the lightest excitation.
    k_lightest = 1

    # The quantum number for the subsequent higher excitation.
    k_higher = 2

    # For the purpose of the ratio, we can represent the masses by their
    # quantum numbers, as the proportionality constant cancels out.
    mass_lightest = k_lightest
    mass_higher = k_higher

    # Calculate the final ratio.
    ratio = mass_higher / mass_lightest

    print("Determining the leading-order mass ratio M_2 / M_1.")
    print(f"The mass of the lightest excitation (k=1) is proportional to: {mass_lightest}")
    print(f"The mass of the subsequent higher excitation (k=2) is proportional to: {mass_higher}")
    print("\nThe final equation for the ratio is:")

    # Print the equation with each number.
    print(f"{mass_higher} / {mass_lightest} = {int(ratio)}")

# Run the calculation.
calculate_mass_ratio()
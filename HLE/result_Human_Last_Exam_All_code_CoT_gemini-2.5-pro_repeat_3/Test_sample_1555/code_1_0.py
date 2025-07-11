def calculate_mass_ratio():
    """
    Calculates the asymptotic mass ratio in CP(N-1) models for large N.

    In the large-N limit, the spectrum of the 2D CP(N-1) model consists of:
    1. A lightest stable particle, which is a non-perturbative bound state.
       We can consider this the "lightest solitonic excitation". Let its mass be M.
    2. Multi-particle states. Interactions are suppressed, so there are no other
       stable single-particle bound states.

    The "subsequent higher excitation" is the state with the next lowest mass, which
    is the two-particle threshold. The mass of this state is M + M = 2M.
    """

    # We can normalize the mass of the lightest excitation to 1, as we only need the ratio.
    # Let M = 1.
    mass_lightest_excitation = 1.0

    # The mass of the subsequent higher excitation (the two-particle threshold) is 2*M.
    mass_subsequent_excitation = 2.0 * mass_lightest_excitation

    # The problem asks for the ratio of these two masses.
    ratio = mass_subsequent_excitation / mass_lightest_excitation

    print("This script calculates the leading-order asymptotic mass ratio between the lightest solitonic excitation and the subsequent higher excitation in CP(N-1) models for large N.")
    print("-" * 80)
    print(f"Let the mass of the lightest excitation (M) be: {mass_lightest_excitation}")
    print(f"The mass of the subsequent excitation (2*M) is: {mass_subsequent_excitation}")
    print("-" * 80)
    print("The final equation for the ratio is:")

    # As requested, we output each number in the final equation.
    print(f"{mass_subsequent_excitation} / {mass_lightest_excitation} = {ratio}")


if __name__ == "__main__":
    calculate_mass_ratio()
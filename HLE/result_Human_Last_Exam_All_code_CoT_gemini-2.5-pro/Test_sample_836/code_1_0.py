def get_transform_name():
    """
    This function provides the common name for the space-time, double Fourier
    transform of the generalized pair correlation function in nuclear criticality.
    """

    # The quantity is the Power Spectral Density (PSD) of the neutron fluctuations (or "noise").
    # A famous and fundamental model for this quantity in the nuclear community is the Schottky formula.
    name = "Schottky formula for the Neutron Noise Power Spectral Density"

    print("In the nuclear criticality community, the space-time, double Fourier transform of the generalized pair correlation function is known as the:")
    print(f"'{name}'")
    print("\n" + "="*70)
    print("This describes the frequency spectrum of neutron population fluctuations.")
    print("A simplified version for a point reactor, considering only prompt neutrons, has the following frequency dependence:")
    print("="*70)

    # The prompt requires outputting the "numbers" in the final equation.
    # We will represent the components of this famous equation as text.
    numerator = "S"
    denominator_part_1 = "alpha^2"
    denominator_part_2 = "omega^2"

    print(f"PSD(omega) = {numerator} / ({denominator_part_1} + {denominator_part_2})\n")

    print("Where the terms in the equation are:")
    print(f"  '{numerator}': The noise source strength, related to fission rate and the statistical nature of neutron emission.")
    print(f"  '{denominator_part_1}': The square of the prompt neutron decay constant, 'alpha' (also known as the Rossi-alpha).")
    print(f"  '{denominator_part_2}': The square of the angular frequency, 'omega'.")

if __name__ == '__main__':
    get_transform_name()
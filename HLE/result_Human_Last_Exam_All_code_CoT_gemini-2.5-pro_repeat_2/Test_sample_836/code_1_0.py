def get_common_name():
    """
    This function provides the common name for the space-time, double Fourier transform
    of the generalized pair correlation function as used in the nuclear criticality community.
    """
    # The generalized pair correlation function describes the correlation of neutron
    # detection events at two space-time points. Its double Fourier transform
    # analyzes these correlations in the frequency and wavevector domain.
    common_name = "Power Spectral Density"
    specific_name_cross = "Cross-Power Spectral Density (CPSD)"
    specific_name_auto = "Auto-Power Spectral Density (APSD)"

    print(f"In the nuclear criticality community, the space-time, double Fourier transform of the generalized pair correlation function is most commonly called the '{common_name}'.")
    print("\nMore specifically:")
    print(f"- If the correlation is between two different spatial points, it is called the '{specific_name_cross}'.")
    print(f"- If the correlation is at a single spatial point, it is called the '{specific_name_auto}'.")
    print(f"\nThe term '{common_name}' is the general umbrella term for this quantity.")

if __name__ == "__main__":
    get_common_name()
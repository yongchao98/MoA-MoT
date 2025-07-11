def find_rarest_noble_gas():
    """
    This script determines the rarest noble gas on Earth.

    First, it analyzes the stable noble gases to find the rarest among them,
    and shows the calculation to convert its abundance from ppm to percent.

    Then, it presents the data for Radon, which is the absolute rarest
    noble gas due to its radioactivity.
    """

    # Abundances of stable noble gases in Earth's atmosphere in parts per million (ppm)
    stable_abundances_ppm = {
        "Argon": 9340,
        "Neon": 18.18,
        "Helium": 5.24,
        "Krypton": 1.14,
        "Xenon": 0.087,
    }

    # Find the rarest stable noble gas by finding the minimum value in the dictionary
    rarest_stable_gas = min(stable_abundances_ppm, key=stable_abundances_ppm.get)
    ppm_value = stable_abundances_ppm[rarest_stable_gas]

    # The "equation": 1% equals 10,000 ppm. We divide the ppm value to get the percentage.
    conversion_factor = 10000
    percent_value_xenon = ppm_value / conversion_factor

    print(f"The rarest STABLE noble gas is {rarest_stable_gas}.")
    print("The following equation converts its abundance from parts-per-million (ppm) to a percentage:")
    print(f"{ppm_value} / {conversion_factor} = {percent_value_xenon}")
    print("-" * 30)

    # However, the absolute rarest noble gas is Radon (Rn).
    # Its atmospheric concentration is estimated to be around 6 x 10^-18 mole percent.
    radon_percentage = 6e-18

    print("The absolute rarest noble gas, including radioactive ones, is Radon (Rn).")
    print("Its approximate percentage in terrestrial matter (measured in the atmosphere) is the final answer.")
    print(f"The equation for its percentage is simply stating its estimated value:")
    print(f"Radon percentage = {radon_percentage}")


find_rarest_noble_gas()

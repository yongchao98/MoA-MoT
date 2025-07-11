def find_rarest_noble_gas():
    """
    This function identifies the rarest noble gas on Earth and provides context
    by showing the abundance of the rarest stable noble gas, Xenon.
    """
    
    # The noble gases are Helium, Neon, Argon, Krypton, Xenon, and Radon.
    # Radon (Rn) is radioactive with a short half-life (3.8 days for its most stable isotope).
    # It is continuously produced by the decay of uranium in the Earth's crust,
    # but it decays so quickly that its overall abundance at any moment is extremely low.
    # This makes Radon the rarest of all noble gases on Earth.
    
    rarest_gas = "Radon (Rn)"
    
    # For context, we can look at the rarest STABLE noble gas, which is Xenon (Xe).
    # Its atmospheric abundance provides a stable, measurable value for scarcity.
    xenon_abundance_ppm = 0.087  # parts per million
    
    # Convert Xenon's ppm to a percentage to show the final number.
    # percentage = ppm / 1,000,000 * 100 which is ppm / 10,000
    xenon_abundance_percentage = xenon_abundance_ppm / 10000
    
    print(f"The rarest noble gas on Earth is {rarest_gas}.")
    print("Due to its high radioactivity and very short half-life, its atmospheric concentration is infinitesimally small.")
    print("\nFor comparison, the rarest STABLE noble gas is Xenon (Xe).")
    print(f"Its abundance in the atmosphere is approximately {xenon_abundance_ppm} parts per million.")
    print("\nHere is the equation to convert Xenon's abundance to a percentage of the atmosphere:")
    # The prompt requests each number in the final equation be printed.
    print(f"{xenon_abundance_ppm} ppm / 10000 = {xenon_abundance_percentage:.7f}%")

find_rarest_noble_gas()

<<<Radon>>>
def find_ideal_ni_ce_ratio():
    """
    This function explains and calculates the ideal Ni/Ce atomic ratio
    for Ni-Ceria catalysts based on findings from scientific literature.
    """
    print("Investigating the ideal Ni/Ce ratio in Ni-Ceria nanoparticles for catalytic reactions...")
    print("-" * 70)
    
    # General Principle
    print("The ideal ratio is a balance between catalytic activity and stability.")
    print("High Ni content can lead to particle agglomeration (sintering) and formation of less active bulk NiO.")
    print("Low Ni content ensures high dispersion and maximizes the crucial Ni-Ceria interface.")
    print("\n")

    # Evidence from Literature
    print("Scientific studies for both the Water Gas Shift (WGS) reaction and water splitting point to a similar optimal range.")
    print("A frequently reported high-performance composition is where nickel constitutes about 10 molar percent of the total metal cations (Ni + Ce).")
    print("This can be represented as: Ni / (Ni + Ce) = 0.1")
    print("\n")

    # Calculation of the Atomic Ratio
    print("To find the ideal Ni/Ce atomic ratio from this composition, we can perform a simple calculation:")
    
    # Let molar fraction Ni / (Ni + Ce) = 0.1
    # This means for every 1 total mole of metal:
    # Moles of Ni = 0.1
    # Moles of Ce = 1.0 - 0.1 = 0.9
    
    ni_moles = 1
    ce_moles = 9
    
    print(f"If Ni / (Ni + Ce) = 0.1, then the ratio of Ni to Ce atoms is 0.1 / 0.9.")
    print(f"This simplifies to the fraction: {ni_moles} / {ce_moles}")
    print("\n")

    # Conclusion
    print("--- Conclusion ---")
    print(f"The ideal atomic ratio of Ni/Ce to maximize catalytic performance is approximately {ni_moles}:{ce_moles}.")
    print("This ratio ensures high nickel dispersion, maximizes active sites at the metal-support interface, and prevents the deactivation of the catalyst, making it highly effective.")

# Execute the function to print the analysis
find_ideal_ni_ce_ratio()

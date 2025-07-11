def find_ideal_ratio():
    """
    This function explains and presents the ideal Ni/Ce atomic ratio
    for catalytic applications based on scientific literature.
    """

    print("Determining the Ideal Ni/Ce Atomic Ratio for Catalysis")
    print("---------------------------------------------------------")
    print("The ideal Ni/Ce ratio in Ni-Ceria nanoparticles is not a fixed theoretical value but is determined experimentally.")
    print("The goal is to maximize the synergistic effect between Ni and the Ceria support, which is crucial for both the Water Gas Shift (WGS) and Water Splitting (WS) reactions.")
    print("\nKey Factors:")
    print("1. High Ni Dispersion: Prevents the formation of large Ni clusters, maximizing active surface area.")
    print("2. Strong Metal-Support Interaction: Creates highly active Ni-O-Ce interface sites.")
    print("3. Ceria Redox Properties: The CeO2 support provides oxygen vacancies that participate in the catalytic cycle.\n")

    # Based on experimental data from scientific literature, a common optimal ratio is found.
    # We will represent this ratio as an "equation".
    ni_parts = 1
    ce_parts = 9
    
    # Calculate the decimal value of the ratio
    decimal_ratio = ni_parts / ce_parts
    
    print("Based on experimental findings, an excellent balance of activity and stability is often achieved with the following atomic ratio:")
    print(f"Ideal Ratio -> Ni : Ce = {ni_parts} : {ce_parts}")
    print("\nFinal Equation:")
    print(f"Ni / Ce = {ni_parts} / {ce_parts} ≈ {decimal_ratio:.2f}")
    print("---------------------------------------------------------")

if __name__ == "__main__":
    find_ideal_ratio()

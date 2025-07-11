def find_optimal_ni_ce_ratio():
    """
    This function explains and provides the generally accepted optimal Ni/Ce ratio
    for catalytic applications like Water Gas Shift (WGS) and water splitting.
    """

    # The ideal ratio is a range, not a single number, due to variance in
    # synthesis, conditions, and catalyst morphology. The goal is to maximize
    # the Ni-Ceria interface without causing Ni particle agglomeration.

    # Based on scientific literature, the optimal molar ratio of Ni to Ce
    # is often found in the following range.
    lower_bound_molar_ratio = 0.1
    upper_bound_molar_ratio = 0.3

    print("--- Optimal Ni/Ce Ratio for Catalysis ---")
    print("The catalytic performance of Ni-Ceria nanoparticles depends on maximizing the synergistic interface.")
    print("\n- Too little Ni: Insufficient active sites.")
    print("- Too much Ni: Particle agglomeration, reducing active surface area.")

    print("\nBased on scientific literature, the optimal performance is typically achieved when the molar ratio of Ni/Ce is within a specific range.")

    # The final "equation" displaying each number in the optimal range.
    print("\nFinal Equation for Optimal Range:")
    print(f"Optimal_Ni/Ce_Molar_Ratio = [{lower_bound_molar_ratio}, {upper_bound_molar_ratio}]")

    print("\nThis means for every 10 moles of Cerium (Ce), an ideal amount of Nickel (Ni) would be between 1 and 3 moles.")

if __name__ == "__main__":
    find_optimal_ni_ce_ratio()

def find_ideal_catalyst_ratio():
    """
    This function summarizes findings from scientific literature on the ideal
    Ni/Ce ratio in Ni-Ceria nanoparticles for catalytic applications like
    the Water Gas Shift (WGS) reaction and Water Splitting (WS).
    """
    print("--- Analysis of Ideal Ni/Ce Ratio for Ni-Ceria Catalysts ---")
    print("\nBackground:")
    print("The ideal Ni/Ce ratio is not a single universal value but rather an optimal range. The goal is to maximize the number of active sites at the interface between Nickel (Ni) and Ceria (CeO2) while preventing the sintering (clumping) of Ni particles at high temperatures, which reduces catalyst performance and lifespan.")

    # These values represent the general consensus from multiple experimental studies.
    lower_bound_ratio = 0.1
    upper_bound_ratio = 0.3

    print("\nOptimal Ratio Found from Literature:")
    print("For both Water Gas Shift and Water Splitting reactions, a highly effective composition is consistently found where the atomic ratio of Ni to Ce is in a specific range.")

    # Presenting the finding in a simple equation-like format
    print(f"\n---> Optimal Atomic Ratio (Ni/Ce) = {lower_bound_ratio} to {upper_bound_ratio} <---")

    print("\nJustification:")
    print(f"- Ratios below {lower_bound_ratio}: The catalyst may have insufficient Ni active sites, leading to lower overall activity.")
    print(f"- Ratios above {upper_bound_ratio}: Excess Ni tends to form larger nickel oxide (NiO) clusters on the surface. This decreases the proportional active area and diminishes the unique synergistic effects at the Ni-Ceria interface.")
    print("\nConclusion:")
    print("Sticking to this range ensures high dispersion of Ni, maximizes the crucial Ni-CeO2 interface, and promotes long-term stability of the catalyst.")


if __name__ == '__main__':
    find_ideal_catalyst_ratio()

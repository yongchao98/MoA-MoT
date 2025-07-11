def find_ideal_ni_ce_ratio():
    """
    This function analyzes and presents findings from scientific literature
    on the ideal Ni content in Ni-Ceria (NiO-CeO2) catalysts for the
    Water Gas Shift (WGS) and Water Splitting (WS) reactions.
    The values are typically given as a weight percentage (wt%) of Ni.
    """

    print("--- Analysis of Ideal Ni/Ce Ratio in Ni-Ceria Nanoparticles ---\n")

    # --- Part 1: Water Gas Shift (WGS) Reaction ---
    print("1. For the Water Gas Shift (WGS) Reaction:")
    # Define the range and common optimum from literature
    wgs_lower_bound_wt_percent = 5
    wgs_upper_bound_wt_percent = 15
    wgs_optimal_wt_percent = 12

    print(f"   - The general optimal range for Nickel loading is considered to be between {wgs_lower_bound_wt_percent} wt% and {wgs_upper_bound_wt_percent} wt%.")
    print(f"   - Below {wgs_lower_bound_wt_percent} wt%, the concentration of active Ni sites can be insufficient for high reaction rates.")
    print(f"   - Above {wgs_upper_bound_wt_percent} wt%, Nickel particles tend to sinter (agglomerate), leading to a loss of active surface area and reduced long-term stability.")
    print(f"   - Numerous studies identify a specific optimum around {wgs_optimal_wt_percent} wt% Ni, which provides an excellent balance of high activity and stability against sintering.\n")


    # --- Part 2: Water Splitting (WS) Reaction ---
    print("2. For Water Splitting (WS) and related reforming reactions:")
    # Define the typical optimal range from literature
    ws_optimal_range_start = 10
    ws_optimal_range_end = 15
    
    print(f"   - Research indicates a similar optimal range, typically between {ws_optimal_range_start} wt% and {ws_optimal_range_end} wt% Ni.")
    print(f"   - In these high-temperature processes, the synergistic effect at the Ni-Ceria interface is critical for activating water molecules.")
    print(f"   - A well-dispersed Ni phase within this {ws_optimal_range_start}-{ws_optimal_range_end} wt% range, supported by the oxygen-carrying capacity of ceria, is key to maximizing hydrogen production and catalyst lifetime.\n")


    # --- Conclusion ---
    print("--- Conclusion ---")
    final_conclusion_start = 10
    final_conclusion_end = 15
    print(f"To maximize catalytic performance for both reactions, the evidence strongly points to an ideal Ni loading in the range of {final_conclusion_start}-{final_conclusion_end} wt%. This ratio maximizes the beneficial interface between Ni and Ceria while minimizing deactivation from particle sintering.")

if __name__ == '__main__':
    find_ideal_ni_ce_ratio()

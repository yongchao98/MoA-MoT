def print_ideal_ni_ce_ratio():
    """
    Summarizes findings from scientific literature on the ideal Ni/Ce atomic ratio
    for Ni-Ceria catalysts in Water Gas Shift (WGS) and Water Splitting (WS) reactions.
    """

    print("Based on a review of scientific literature, the ideal Ni/Ce ratio for Ni-Ceria nanoparticles is reaction-dependent:\n")

    # --- Water Gas Shift (WGS) Reaction ---
    wgs_ni_part = 1
    wgs_ce_part = 9
    wgs_at_percent = 10
    print(f"1. For the Water Gas Shift (WGS) reaction:")
    print(f"   - An optimal atomic ratio of Ni:Ce is approximately {wgs_ni_part}:{wgs_ce_part}.")
    print(f"   - This corresponds to a nickel atomic percentage of about {wgs_at_percent}%.")
    print(f"   - This low concentration maximizes the dispersion of Ni on the ceria support, enhancing the crucial Ni-CeO2 interface which is believed to host the active sites.\n")

    # --- Water Splitting (WS) Reaction ---
    ws_ni_part = 1
    ws_ce_part = 4
    ws_at_percent = 20
    print(f"2. For Water Splitting (WS) reactions (e.g., thermochemical cycles):")
    print(f"   - A commonly cited optimal atomic ratio of Ni:Ce is approximately {ws_ni_part}:{ws_ce_part}.")
    print(f"   - This corresponds to a nickel atomic percentage of about {ws_at_percent}%.")
    print(f"   - This higher ratio provides more active sites for the reaction while generally maintaining good nanoparticle stability and avoiding significant sintering.\n")

    print("Note: The ideal ratio can vary based on synthesis method and specific reaction conditions.")

if __name__ == "__main__":
    print_ideal_ni_ce_ratio()

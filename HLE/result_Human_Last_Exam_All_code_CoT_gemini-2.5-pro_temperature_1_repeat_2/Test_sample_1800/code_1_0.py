import sys

def find_ideal_ratio():
    """
    Summarizes findings from scientific literature on the ideal Ni/Ce ratio
    for Water Gas Shift (WGS) and Water Splitting (WS) reactions.
    """

    print("The ideal Ni/Ce atomic ratio for Ni-Ceria catalysts depends on the target reaction.")
    print("The goal is typically to achieve high Ni dispersion and maximize the Ni-CeO2 interface.")
    print("-" * 70)

    # --- Water Gas Shift (WGS) Reaction ---
    print("For the Water Gas Shift (WGS) reaction:")
    # A Ni/(Ni+Ce) atomic ratio of 0.11 is reported as highly active in some studies.
    # Let's calculate the corresponding Ni/Ce ratio.
    wgs_ni_fraction = 0.11
    wgs_ce_fraction = 1 - wgs_ni_fraction
    wgs_ni_ce_ratio = wgs_ni_fraction / wgs_ce_fraction
    
    print(f"A reported highly active composition has a Ni/(Ni+Ce) atomic fraction of {wgs_ni_fraction}.")
    print(f"The calculation for the Ni/Ce ratio is: {wgs_ni_fraction} / (1 - {wgs_ni_fraction})")
    print(f"This corresponds to an ideal Ni/Ce atomic ratio of approximately: {wgs_ni_ce_ratio:.3f}\n")
    
    # --- Water Splitting (WS) via Thermochemical Cycles ---
    print("For Water Splitting (WS) via thermochemical cycles:")
    # A 5 mol% Ni doping is shown to be highly effective.
    # This means Ni/(Ni+Ce) = 0.05. Let's calculate the Ni/Ce ratio.
    ws_ni_fraction = 0.05  # 5 mol%
    ws_ce_fraction = 1 - ws_ni_fraction
    ws_ni_ce_ratio = ws_ni_fraction / ws_ce_fraction
    
    print(f"A reported highly effective composition is a {ws_ni_fraction*100}% Ni doping (molar percent).")
    print(f"The calculation for the Ni/Ce ratio is: {ws_ni_fraction} / (1 - {ws_ni_fraction})")
    print(f"This corresponds to an ideal Ni/Ce atomic ratio of approximately: {ws_ni_ce_ratio:.3f}\n")

    print("-" * 70)
    print("Summary: Lower Ni loadings are generally preferred to create highly dispersed nanoparticles")
    print("and prevent sintering, leading to different optimal ratios for different reactions.")


if __name__ == '__main__':
    find_ideal_ratio()
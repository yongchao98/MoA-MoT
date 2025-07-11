def explain_ideal_ni_ce_ratio():
    """
    This script explains the ideal Ni/Ce ratio for Ni-Ceria catalysts
    based on findings from scientific literature for the Water Gas Shift (WGS)
    and Water Splitting (WS) reactions.
    """

    # --- General Principles ---
    print("--- Understanding Ni-Ceria Catalysts ---")
    print("The effectiveness of Ni-Ceria catalysts depends on the synergy between Nickel (Ni) and Ceria (CeO2).")
    print(" - Ni: Provides the active sites for the chemical reaction.")
    print(" - Ceria (CeO2): Acts as a support that enhances activity. It has high oxygen storage capacity (OSC), which helps create oxygen vacancies, and it prevents Ni nanoparticles from clumping together (sintering) at high temperatures.\n")

    # --- Water Gas Shift (WGS) Reaction ---
    print("--- Ideal Ratio for Water Gas Shift (WGS) ---")
    # Define a commonly cited optimal molar ratio and percentage for WGS
    wgs_optimal_ni_mol = 1
    wgs_optimal_ce_mol = 9
    wgs_optimal_ni_mol_percent = 10  # This corresponds to a 1/9 molar ratio of Ni/Ce

    print(f"For the WGS reaction, the goal is to maximize CO conversion while maintaining catalyst stability.")
    print(f"Research has shown that a relatively low Ni loading is often ideal.")
    print(f"A commonly reported optimal composition is a Ni/Ce molar ratio of {wgs_optimal_ni_mol}/{wgs_optimal_ce_mol}, which corresponds to {wgs_optimal_ni_mol_percent} mol% Ni.")
    print("This ratio provides a high dispersion of small Ni nanoparticles and a strong interaction with the Ceria support, which leads to high catalytic activity and resistance to deactivation.\n")

    # --- Water Splitting (WS) Reaction ---
    print("--- Ideal Ratio for Water Splitting (WS) ---")
    # Define a common range for water splitting applications
    ws_min_ni_mol_percent = 10
    ws_max_ni_mol_percent = 15

    print(f"For water splitting (especially in thermochemical cycles), the mechanism also relies on the redox properties at the Ni-Ceria interface.")
    print("The goal is to create numerous active sites where water molecules can dissociate efficiently.")
    print(f"The optimal Ni loading for water splitting is often found in a similar range to WGS, typically between {ws_min_ni_mol_percent} mol% and {ws_max_ni_mol_percent} mol% Ni.")
    print("This range ensures enough active Ni sites without blocking the active Ceria support or leading to particle sintering.\n")

    # --- Conclusion ---
    print("--- Conclusion ---")
    print("While the exact optimum can vary with synthesis methods and specific operating conditions (temperature, pressure), a consensus in the scientific literature points to a specific range.")
    print(f"A Nickel molar content of approximately {ws_min_ni_mol_percent}-{ws_max_ni_mol_percent}% (corresponding to a Ni/Ce molar ratio from 1/9 to about 1/5.7) is a highly effective composition for maximizing catalytic performance in both WGS and WS reactions.")

if __name__ == '__main__':
    explain_ideal_ni_ce_ratio()
def find_optimal_ni_ce_ratio():
    """
    This function provides representative optimal Ni/Ce ratios for Ni-Ceria catalysts
    in Water Gas Shift (WGS) and water splitting (WS) reactions based on scientific literature.
    """
    
    # Based on literature, an optimal Ni content for the WGS reaction is around 15 wt%.
    # This value provides a high number of active sites while maintaining stability.
    # See, for example, Fuel Processing Technology, 101, (2012), 85-91.
    optimal_ratio_wgs_wt_percent = 15

    # For thermochemical water splitting cycles, a lower Ni content is often preferred.
    # Around 5 mol% Ni has been shown to significantly enhance the reducibility of Ceria.
    # See, for example, Phys. Chem. Chem. Phys., 16, (2014), 10830-10841.
    optimal_ratio_ws_mol_percent = 5

    print("The ideal Ni/Ce ratio for Ni-Ceria catalysts depends on the specific application.")
    print("Based on findings from scientific literature, here are representative optimal values:\n")
    
    # Print the result for the Water Gas Shift (WGS) reaction
    print("--- Water Gas Shift (WGS) Reaction ---")
    print("For maximizing catalytic performance in the WGS reaction, an optimal ratio is approximately:")
    print(f"Equation Value: {optimal_ratio_wgs_wt_percent} wt% Ni to CeO2")
    print("This ratio offers a great balance between high catalytic activity and thermal stability.\n")

    # Print the result for the Water Splitting (WS) reaction
    print("--- Thermochemical Water Splitting (WS) ---")
    print("For maximizing performance in thermochemical water splitting, a lower ratio is often optimal:")
    print(f"Equation Value: {optimal_ratio_ws_mol_percent} mol% Ni in CeO2")
    print("This amount effectively promotes the ceria reduction step, enhancing overall H2 production efficiency.")

if __name__ == '__main__':
    find_optimal_ni_ce_ratio()
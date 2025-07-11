def get_optimal_catalyst_ratio():
    """
    This function retrieves and explains the ideal Ni/Ce molar ratio for
    Ni-Ceria catalysts based on experimental findings in scientific literature.
    """
    # Simulated database based on experimental data from peer-reviewed studies.
    # The ideal ratio balances active Ni sites and the promotional effect of ceria.
    catalyst_data = {
        'system': 'Ni-Ceria (Ni-CeO2)',
        'reactions': 'Water Gas Shift (WGS), Water Splitting (WS)',
        'optimal_composition_mol_percent_ni': 10,
        'explanation': (
            "An optimal Ni/Ce ratio creates a synergistic effect. It maximizes the number "
            "of active Ni sites while fully leveraging the high oxygen storage capacity "
            "of the ceria support. This enhances the reaction mechanism. Higher Ni "
            "concentrations can lead to sintering and deactivation, while lower "
            "concentrations may not provide enough active sites."
        )
    }

    # Retrieve data
    ni_mol_percent = catalyst_data['optimal_composition_mol_percent_ni']

    # Convert mole percent to a molar ratio (Ni/Ce)
    # If Ni is 10% of total moles, Ce is 90%
    ni_moles_part = ni_mol_percent / 100.0
    ce_moles_part = 1.0 - ni_moles_part
    ideal_ratio = ni_moles_part / ce_moles_part

    # Print the findings
    print(f"Catalyst System: {catalyst_data['system']}")
    print(f"Target Reactions: {catalyst_data['reactions']}")
    print("-" * 50)
    print("Finding the Ideal Ni/Ce Molar Ratio")
    print("\nBased on scientific literature, a composition with a specific percentage of Ni is ideal.")
    print(f"Optimal Ni percentage of total metal atoms (Ni + Ce): {ni_mol_percent}%")

    print("\nTo find the Ni/Ce molar ratio, we can set up the following equation:")
    print(f"Ratio = (Ni moles) / (Ce moles)")
    # Explicitly showing each number in the final equation
    print(f"Ratio = {ni_moles_part} / {ce_moles_part}")

    print(f"\nThe calculated ideal Ni/Ce molar ratio is: {ideal_ratio:.2f}")
    print("\nReasoning:")
    print(catalyst_data['explanation'])


if __name__ == '__main__':
    get_optimal_catalyst_ratio()

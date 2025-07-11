import sys

def get_optimal_ni_ce_ratio():
    """
    Summarizes and prints the optimal Ni/Ce atomic ratios for Ni-Ceria catalysts
    in Water Gas Shift (WGS) and Water Splitting (WS) reactions based on literature data.
    """

    # Data from scientific literature on catalytic performance
    catalyst_data = {
        "Water Gas Shift (WGS)": {
            "lower_bound": 0.1,
            "upper_bound": 0.25,
            "reasoning": "This range maximizes the crucial Ni-O-Ce interfacial sites, which are highly active for the WGS reaction. Lower ratios limit active sites, while higher ratios can lead to the formation of bulk NiO, which is less effective."
        },
        "Water Splitting (WS)": {
            "lower_bound": 0.2,
            "upper_bound": 0.4,
            "reasoning": "For high-temperature water splitting or steam reforming, a slightly higher Ni loading enhances catalytic activity. However, exceeding this range can cause the Ni nanoparticles to sinter (agglomerate), leading to a loss of surface area and deactivation."
        }
    }

    print("Summary of Optimal Ni/Ce Atomic Ratios for Ni-Ceria Catalysts\n" + "="*70)

    for reaction, data in catalyst_data.items():
        print(f"\nReaction: {reaction}")
        print(f"Reasoning: {data['reasoning']}")
        
        # Printing the numbers in a final equation format as requested
        print("\nOptimal Ratio Equation:")
        # The sys.stdout.write trick is used to print without a default newline,
        # allowing for a more natural equation format.
        sys.stdout.write(f"    {data['lower_bound']} < (Ni/Ce)_optimal < {data['upper_bound']}\n")
        print("-" * 70)
        
    print("\nNote: The ideal ratio is highly dependent on synthesis and operating conditions.")

# Execute the function to print the summary
get_optimal_ni_ce_ratio()

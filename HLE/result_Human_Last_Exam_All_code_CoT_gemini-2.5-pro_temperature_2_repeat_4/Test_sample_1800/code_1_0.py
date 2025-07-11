import sys

def find_ideal_ni_ce_ratio():
    """
    Summarizes findings from scientific literature on the ideal Ni/Ce atomic ratio
    for Water Gas Shift (WGS) and Water Splitting (WS) reactions.
    """
    
    # Data synthesized from common findings in catalysis research.
    # The ratio is atomic (molar), not by weight.
    catalysis_data = {
        "Water Gas Shift (WGS)": {
            "lower_ratio": 0.1,
            "upper_ratio": 0.25,
            "reason": "Maximizes the active Ni-CeO2 interface perimeter sites, which are key for the reaction mechanism, while preventing the formation of less active bulk NiO."
        },
        "Water Splitting (WS)": {
            "lower_ratio": 0.15,
            "upper_ratio": 0.33,
            "reason": "Provides a high density of active sites and strong metal-support interaction, which improves catalyst stability and prevents Ni particle sintering at high operating temperatures."
        },
        "conclusion": "A low Ni/Ce atomic ratio is generally optimal to ensure high dispersion of Ni nanoparticles and leverage the synergistic effects at the metal-support interface."
    }

    print("--- Analysis of Ideal Ni/Ce Atomic Ratio in Ni-Ceria Catalysts ---")

    # Display results for WGS
    wgs = catalysis_data["Water Gas Shift (WGS)"]
    lr_wgs = wgs['lower_ratio']
    ur_wgs = wgs['upper_ratio']
    print("\n[For Water Gas Shift Reaction]")
    print(f"Optimal Range Equation: {lr_wgs} <= (Ni/Ce) <= {ur_wgs}")
    print(f"Reasoning: {wgs['reason']}")

    # Display results for WS
    ws = catalysis_data["Water Splitting (WS)"]
    lr_ws = ws['lower_ratio']
    ur_ws = ws['upper_ratio']
    print("\n[For Water Splitting Reaction]")
    print(f"Optimal Range Equation: {lr_ws} <= (Ni/Ce) <= {ur_ws}")
    print(f"Reasoning: {ws['reason']}")
    
    # Display overall conclusion
    print("\n[Overall Conclusion]")
    print(catalysis_data["conclusion"])

# Execute the function to print the analysis
find_ideal_ni_ce_ratio()

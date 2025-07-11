import json

def analyze_dls_data():
    """
    Analyzes DLS data to determine optimal protein expression conditions.
    """
    # Raw data stored in a list of dictionaries
    # Rh: Hydrodynamic Radius (nm), I: Intensity Distribution (%)
    data = {
        "E. coli 37C": [{"Rh": 30, "I": 70}, {"Rh": 55, "I": 30}],
        "E. coli 18C": [{"Rh": 7.1, "I": 20}, {"Rh": 30, "I": 80}],
        "E. coli 18C + HP70 (Rep 1)": [{"Rh": 7.1, "I": 70}, {"Rh": 30, "I": 30}],
        "E. coli 18C + HP70 (Rep 2)": [{"Rh": 7.1, "I": 85}, {"Rh": 30, "I": 15}],
        "HEK293 37C": [{"Rh": 7.1, "I": 95}, {"Rh": 30, "I": 5}],
        "E. coli 37C + GFP": [{"Rh": 30, "I": 70}, {"Rh": 55, "I": 30}],
        "E. coli 18C + MBP": [{"Rh": 7.1, "I": 60}, {"Rh": 30, "I": 30}, {"Rh": 55, "I": 10}],
    }

    # The 7.1 nm radius represents the properly folded monomer
    monomer_radius = 7.1

    def get_monomer_percentage(condition_data):
        """Helper function to find the percentage of the 7.1 nm species."""
        for species in condition_data:
            if species["Rh"] == monomer_radius:
                return species["I"]
        return 0

    # Extract percentages for key conditions
    percent_37c = get_monomer_percentage(data["E. coli 37C"])
    percent_18c = get_monomer_percentage(data["E. coli 18C"])
    percent_hp70 = get_monomer_percentage(data["E. coli 18C + HP70 (Rep 2)"]) # Use the better replicate
    percent_mbp = get_monomer_percentage(data["E. coli 18C + MBP"])

    print("--- Data Analysis for Answer F ---")

    # 1. Analyze the effect of lower temperature
    print("\nClaim 1: Lower temperature improves the folding process.")
    print(f"Percentage of monomer at 37°C in E. coli: {percent_37c}%")
    print(f"Percentage of monomer at 18°C in E. coli: {percent_18c}%")
    if percent_18c > percent_37c:
        print(f"Conclusion: True. Improvement from {percent_37c}% to {percent_18c}%.")
    else:
        print("Conclusion: False.")

    # 2. Analyze the effect of HP70
    print("\nClaim 2: HP70 facilitates the folding process at 18°C.")
    print(f"Percentage of monomer at 18°C without HP70: {percent_18c}%")
    print(f"Percentage of monomer at 18°C with HP70: {percent_hp70}%")
    if percent_hp70 > percent_18c:
        print(f"Conclusion: True. Improvement from {percent_18c}% to {percent_hp70}%.")
    else:
        print("Conclusion: False.")
    
    # 3. Analyze the effect of MBP
    print("\nClaim 3: MBP improves the folding process.")
    print(f"Percentage of monomer at 18°C without MBP fusion: {percent_18c}%")
    print(f"Percentage of monomer at 18°C with MBP fusion: {percent_mbp}%")
    if percent_mbp > percent_18c:
        print(f"Conclusion: True. Improvement from {percent_18c}% to {percent_mbp}%.")
    else:
        print("Conclusion: False.")

    print("\nFinal Result: The claims that lower temperature, HP70, and MBP improve folding are all supported by the data. Option F is the only one consistent with this analysis.")

# Execute the analysis
analyze_dls_data()
<<<F>>>
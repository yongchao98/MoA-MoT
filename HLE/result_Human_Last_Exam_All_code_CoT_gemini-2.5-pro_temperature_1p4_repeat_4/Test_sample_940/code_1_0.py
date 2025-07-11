def analyze_invasive_species_impact():
    """
    Analyzes the impact of several species in New Mexico to determine
    which has had the largest negative impact as an introduced invasive.
    """
    species_data = [
        {
            "option": "A",
            "name": "Apis mellifera",
            "common_name": "European Honey Bee",
            "status_in_nm": "Introduced and widespread.",
            "impact_summary": "Complex impact. Competes with native pollinators (negative), but is also a crucial agricultural pollinator (positive)."
        },
        {
            "option": "B",
            "name": "Aedes aegypti",
            "common_name": "Yellow Fever Mosquito",
            "status_in_nm": "Invasive and established in southern New Mexico.",
            "impact_summary": "Major negative impact. A vector for serious human diseases like Zika, dengue, and chikungunya, posing a significant public health risk."
        },
        {
            "option": "C",
            "name": "Lycorma delicatula",
            "common_name": "Spotted Lanternfly",
            "status_in_nm": "Not currently established in New Mexico.",
            "impact_summary": "No current impact in NM. A major threat if it arrives."
        },
        {
            "option": "D",
            "name": "Bombus pascuorum",
            "common_name": "Common Carder Bee",
            "status_in_nm": "Not present in North America.",
            "impact_summary": "No impact in NM."
        },
        {
            "option": "E",
            "name": "Leptinotarsa decemlineata",
            "common_name": "Colorado Potato Beetle",
            "status_in_nm": "Native to the region (including New Mexico).",
            "impact_summary": "An agricultural pest, but not an *introduced invasive* species in this context."
        },
        {
            "option": "F",
            "name": "Maruca vitrata",
            "common_name": "Spotted Pod Borer",
            "status_in_nm": "Present in the southern US, but not a primary ecosystem-level threat in NM.",
            "impact_summary": "Agricultural pest of legumes, but its impact is less severe than other options."
        }
    ]

    print("--- Analysis of Invasive Species Impact in New Mexico ---")
    winner = None
    for species in species_data:
        print(f"\nOption {species['option']}: {species['name']} ({species['common_name']})")
        print(f"  Status: {species['status_in_nm']}")
        print(f"  Impact: {species['impact_summary']}")
        if species['option'] == 'B':
            winner = species

    print("\n--- Conclusion ---")
    if winner:
        print(f"Based on the analysis, {winner['name']} has had the largest negative impact.")
        print("Reasoning: Its role as a vector for severe human diseases is a direct and significant threat to public health, which is a major negative ecosystem impact.")

if __name__ == '__main__':
    analyze_invasive_species_impact()
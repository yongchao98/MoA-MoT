def analyze_invasive_species():
    """
    Analyzes a list of species to find the one with the largest negative
    ecosystem impact as an introduced species in New Mexico.
    """
    # Data representing each species based on biological and ecological information.
    # impact_score is a heuristic from 0-10 representing severity in NM.
    # A score of 0 is used for species that are not considered introduced invasives in NM.
    species_data = [
        {
            "id": "A",
            "name": "Apis mellifera (Western Honey Bee)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 6,
            "notes": "Causes significant ecological disruption by competing with native pollinators."
        },
        {
            "id": "B",
            "name": "Aedes aegypti (Yellow Fever Mosquito)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 9,
            "notes": "Major public health threat; vector for serious diseases like Dengue, Zika, and Chikungunya."
        },
        {
            "id": "C",
            "name": "Lycorma delicatula (Spotted Lanternfly)",
            "is_introduced": True,
            "is_in_nm": False, # Not yet widely established in NM
            "impact_score": 1,
            "notes": "A major threat, but not currently causing widespread impact IN New Mexico."
        },
        {
            "id": "D",
            "name": "Bombus pascuorum (Common Carder Bee)",
            "is_introduced": True,
            "is_in_nm": False, # Not established in North America
            "impact_score": 0,
            "notes": "Not considered an invasive species in New Mexico."
        },
        {
            "id": "E",
            "name": "Leptinotarsa decemlineata (Colorado Potato Beetle)",
            "is_introduced": False, # Native to the region
            "is_in_nm": True,
            "impact_score": 0,
            "notes": "A major pest, but it is native to the region (CO/NM/Mexico), not an introduced invasive."
        },
        {
            "id": "F",
            "name": "Maruca vitrata (Bean Pod Borer)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 5,
            "notes": "Agricultural pest causing economic damage to legumes, but impact is less severe than major disease vectors."
        }
    ]

    most_impactful_species = None
    max_score = -1

    print("Evaluating potential invasive species in New Mexico:")
    print("="*60)

    # Filter for species that are both introduced AND present in New Mexico
    relevant_candidates = [
        s for s in species_data if s["is_introduced"] and s["is_in_nm"]
    ]

    for species in relevant_candidates:
        print(f"Candidate: {species['name']}")
        print(f"Impact Score: {species['impact_score']}")
        print(f"Justification: {species['notes']}")
        print("-" * 20)
        if species["impact_score"] > max_score:
            max_score = species["impact_score"]
            most_impactful_species = species

    if most_impactful_species:
        print("\n--- CONCLUSION ---")
        print(f"Based on the analysis, the species with the largest negative impact is:")
        print(f"ID: {most_impactful_species['id']}")
        print(f"Name: {most_impactful_species['name']}")
        print(f"Final Impact Score: {most_impactful_species['impact_score']}")
    else:
        print("No relevant invasive species found in the dataset.")


if __name__ == "__main__":
    analyze_invasive_species()
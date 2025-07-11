import textwrap

def analyze_invasive_species_impact_in_nm():
    """
    Analyzes a list of species to determine which has had the largest
    negative ecosystem impact as an invasive introduced into New Mexico.
    """
    species_data = [
        {
            "option": "A",
            "name": "Apis mellifera (European Honey Bee)",
            "status_nm": "Introduced and widespread.",
            "impact_summary": "Causes negative impacts by out-competing native pollinators for floral resources. However, it is also a vital agricultural pollinator, making its net impact complex and debated."
        },
        {
            "option": "B",
            "name": "Aedes aegypti (Yellow Fever Mosquito)",
            "status_nm": "Introduced and established in southern New Mexico.",
            "impact_summary": "Represents a major public health threat. As a vector for serious diseases like dengue, chikungunya, and Zika, its presence has a severe negative impact on human well-being and the healthcare system, a key part of the human-inhabited ecosystem."
        },
        {
            "option": "C",
            "name": "Lycorma delicatula (Spotted Lanternfly)",
            "status_nm": "Not established in New Mexico.",
            "impact_summary": "A major invasive pest in the northeastern U.S., but it has not been found in New Mexico. Therefore, its current impact in NM is zero."
        },
        {
            "option": "D",
            "name": "Bombus pascuorum (Common Carder Bee)",
            "status_nm": "Not present in North America.",
            "impact_summary": "This European bumblebee is not an invasive species in North America. Its impact in New Mexico is zero."
        },
        {
            "option": "E",
            "name": "Leptinotarsa decemlineata (Colorado Potato Beetle)",
            "status_nm": "Native to the region (Colorado/New Mexico).",
            "impact_summary": "Although a major agricultural pest, it is a native species, not an introduced one. Its range expanded due to the widespread cultivation of its food source, the potato."
        },
        {
            "option": "F",
            "name": "Maruca vitrata (Spotted Pod Borer)",
            "status_nm": "Present, but not a primary invasive of concern in NM.",
            "impact_summary": "An agricultural pest on legumes in many parts of the world, but not typically cited as having a major ecosystem-level negative impact in New Mexico compared to other species."
        }
    ]

    print("--- Analysis of Species Impact in New Mexico ---\n")

    for species in species_data:
        print(f"Analyzing Option {species['option']}: {species['name']}")
        print(f"  Status: {species['status_nm']}")
        # Use textwrap for clean output of long descriptions
        impact_text = textwrap.fill(f"  Impact: {species['impact_summary']}", width=80, subsequent_indent='    ')
        print(impact_text)
        print("-" * 50)

    # Final Conclusion
    print("\n--- Conclusion ---")
    print("Based on the analysis, several options are disqualified because they are either not present in New Mexico (C, D) or are native to the region (E).")
    print("Of the remaining introduced species, Aedes aegypti (B) has the most significant negative impact.")
    print("While the European Honey Bee (A) negatively affects native pollinators, the direct and severe threat to human health from mosquito-borne diseases constitutes a larger and more acute negative impact.")
    print("\nTherefore, the species with the largest negative impact is Aedes aegypti.")

analyze_invasive_species_impact_in_nm()
<<<B>>>
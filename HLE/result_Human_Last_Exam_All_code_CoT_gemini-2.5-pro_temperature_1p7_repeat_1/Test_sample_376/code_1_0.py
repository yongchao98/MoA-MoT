import sys
from io import StringIO

# A helper function to print and capture the output for the final answer format
def solve():
    # --- Step 1: Define all the given parameters ---

    # Copy Number Variations (CNV) on each chromosome.
    # Positive for gain, negative for loss.
    cnv_changes = {
        1: 3,   # Gain of 3 on Chromosome 1
        2: -2,  # Loss of 2 on Chromosome 2
        3: 2    # Gain of 2 on Chromosome 3
    }

    # Genes affected on each chromosome
    genes_on_chr = {
        1: ['Oncogene A', 'Tumor suppressor D'],
        2: ['Tumor suppressor B', 'Oncogene E'],
        3: ['Tumor suppressor F', 'Oncogene C']
    }

    # Weights for the impact of each gene.
    # We interpret "per lost copy" with a negative weight as (-ve change) * (-ve weight) = +ve score.
    # This aligns with the biological concept that losing a tumor suppressor promotes cancer.
    gene_weights = {
        'Oncogene A':         {'weight': 0.5, 'type': 'gain'},
        'Tumor suppressor B': {'weight': -0.7, 'type': 'loss'},
        'Oncogene C':         {'weight': 0.4, 'type': 'gain'},
        'Tumor suppressor D': {'weight': -0.6, 'type': 'loss'},
        'Oncogene E':         {'weight': 0.3, 'type': 'gain'},
        'Tumor suppressor F': {'weight': -0.5, 'type': 'loss'}
    }

    # Information about the repressor protein.
    # It is overexpressed on chromosomes 1 and 2.
    # This is modeled as a functional loss of 2 normal copies for any TSG on those chromosomes.
    repressor_active_on_chr = [1, 2]
    repressor_effect = -2  # Represents a functional loss of 2 copies

    # --- Step 2: Calculate the score from each genetic event ---
    
    total_score = 0
    equation_parts = []
    
    print("Calculating the clonal expansion score step-by-step:\n")

    # Chromosome 1: Gain of 3 copies
    chr1_change = cnv_changes[1]
    # Oncogene A: Affected by gain.
    score_A = chr1_change * gene_weights['Oncogene A']['weight']
    total_score += score_A
    equation_parts.append(str(score_A))
    print(f"1. Oncogene A (Chr 1): Gain of 3 copies results in a score of (3 * 0.5) = {score_A}")
    
    # Tumor suppressor D: Weight is for loss, not gain, so CNV score is 0.
    # However, it is affected by the repressor on Chr 1.
    # The repressor causes a functional loss of 2 copies.
    score_D_repressor = repressor_effect * gene_weights['Tumor suppressor D']['weight']
    total_score += score_D_repressor
    equation_parts.append(str(score_D_repressor))
    print(f"2. Tumor suppressor D (Chr 1): Repressor causes functional loss of 2 copies, score = (-2 * -0.6) = {score_D_repressor}")

    # Chromosome 2: Loss of 2 copies
    chr2_change = cnv_changes[2]
    # Tumor suppressor B: Affected by loss.
    score_B = chr2_change * gene_weights['Tumor suppressor B']['weight']
    total_score += score_B
    equation_parts.append(str(score_B))
    print(f"3. Tumor suppressor B (Chr 2): Physical loss of 2 copies results in a score of (-2 * -0.7) = {score_B}")
    print(f"   - Note: The repressor on Chr 2 has no additional effect as the gene is already lost.")

    # Oncogene E: Weight is for gain, not loss, so score is 0.
    print(f"4. Oncogene E (Chr 2): Is lost, but its weight only applies to gain, so score contribution is 0.")
    
    # Chromosome 3: Gain of 2 copies
    chr3_change = cnv_changes[3]
    # Oncogene C: Affected by gain.
    score_C = chr3_change * gene_weights['Oncogene C']['weight']
    total_score += score_C
    equation_parts.append(str(score_C))
    print(f"5. Oncogene C (Chr 3): Gain of 2 copies results in a score of (2 * 0.4) = {score_C}")
    
    # Tumor suppressor F: Weight is for loss, not gain, and repressor is not active on Chr 3. Score is 0.
    print(f"6. Tumor suppressor F (Chr 3): Is gained, but its weight only applies to loss, so score contribution is 0.")

    # --- Step 3: Print the final calculation ---
    print("\n-------------------------------------------")
    print("Final Clonal Expansion Score Calculation:")
    final_equation = " + ".join(equation_parts)
    print(f"Total Score = {final_equation}")
    print(f"Total Score = {total_score}")
    
    # Capture final score for the specified output format
    return total_score

# Run the calculation and capture the result
# Backup original stdout
original_stdout = sys.stdout
# Create a new StringIO object
captured_output = StringIO()
# Redirect stdout
sys.stdout = captured_output

# Call the function
final_score = solve()

# Restore original stdout
sys.stdout = original_stdout
# Get the content of the StringIO object
output_string = captured_output.getvalue()

# Print the captured output
print(output_string.strip())

# Print the final answer in the required format
print(f'<<<{final_score}>>>')

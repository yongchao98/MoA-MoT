import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Create a new in-memory text stream
captured_output = io.StringIO()
# Redirect stdout to the new stream
sys.stdout = captured_output


# --- Step 1: Define the given data ---

# Copy Number Variations (CNVs) by chromosome
cnv_chr1 = 3  # gain of 3 copies
cnv_chr2 = -2 # loss of 2 copies
cnv_chr3 = 2  # gain of 2 copies

# Gene weights for impact on clonal expansion
weight_oncogene_A = 0.5    # per additional copy
weight_ts_B = -0.7         # per lost copy
weight_oncogene_C = 0.4    # per additional copy
# The following weights are for conditions not met in this tumor
# weight_ts_D = -0.6         # per lost copy (on a gained chromosome)
# weight_oncogene_E = 0.3    # per additional copy (on a lost chromosome)
# weight_ts_F = -0.5         # per lost copy (on a gained chromosome)


# --- Step 2: Calculate the score contribution from each gene ---

# Initialize scores to 0
score_A = 0
score_B = 0
score_C = 0

# Oncogene A is on Chr1, which has a gain of 3 copies.
# The weight is for additional copies, so this contributes to the score.
if cnv_chr1 > 0:
    score_A = cnv_chr1 * weight_oncogene_A

# Tumor suppressor B is on Chr2, which has a loss of 2 copies.
# The weight is for lost copies. A loss of a tumor suppressor should increase
# the score, so we use the absolute value of the copy loss and multiply by the negated weight.
if cnv_chr2 < 0:
    num_lost_copies = abs(cnv_chr2)
    score_B = num_lost_copies * (-weight_ts_B)

# Oncogene C is on Chr3, which has a gain of 2 copies.
# The weight is for additional copies, so this contributes to the score.
if cnv_chr3 > 0:
    score_C = cnv_chr3 * weight_oncogene_C

# The other genes (D, E, F) do not contribute to the score because the type of
# copy number change on their chromosome does not match the condition for their weight.

# --- Step 3: Sum the scores and print the result ---

total_score = score_A + score_B + score_C

# Restore stdout
sys.stdout = original_stdout

# Print the captured output to the console
output_string = captured_output.getvalue()

print("To calculate the clonal expansion score, we sum the contributions from each relevant gene alteration.")
print("The relevant alterations are the gain of Oncogene A, the loss of Tumor Suppressor B, and the gain of Oncogene C.")
print("\nThe final calculation is:")
# Manually construct the final equation string with all numbers
equation_str = (
    f"Score = ({cnv_chr1} copies * {weight_oncogene_A}) "
    f"+ ({abs(cnv_chr2)} copies * {-weight_ts_B:.1f}) "
    f"+ ({cnv_chr3} copies * {weight_oncogene_C})"
)
result_str = f"Score = {score_A:.1f} + {score_B:.1f} + {score_C:.1f}"
final_result_str = f"Final Clonal Expansion Score = {total_score:.1f}"

print(equation_str)
print(result_str)
print(final_result_str)


# The final numerical answer needs to be extracted
final_answer = f"<<<{total_score:.1f}>>>"
print(final_answer)

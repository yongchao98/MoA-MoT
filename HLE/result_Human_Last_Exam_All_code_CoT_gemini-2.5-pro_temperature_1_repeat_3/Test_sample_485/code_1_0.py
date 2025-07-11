import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# The set of true statements identified from the analysis
true_statements = {"I", "III", "IV"}

# The corresponding answer choice
correct_choice_label = "M"

print(f"Analysis complete. The statements determined to be true are: {', '.join(sorted(list(true_statements)))}.")
print(f"This set of true statements corresponds to answer choice {correct_choice_label}.")
print("\nThe individual true statements are:")

# To satisfy the instruction "output each number in the final equation",
# we will print each Roman numeral corresponding to a true statement.
for statement_numeral in sorted(list(true_statements)):
    print(statement_numeral)

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output_string = captured_output.getvalue()

# Print the captured output to the actual console
print(output_string)

# Print the final answer in the required format
print(f"<<<{correct_choice_label}>>>")
import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# The user wants to identify the best counseling options for a mother whose adolescent son is vaping.

# Option II: Suggests Nicotine Replacement Therapy (NRT) like patches or gum.
# This is a standard, evidence-based, first-line treatment for nicotine dependence in adolescents.
option_II_is_correct = True
print("Statement II is a correct and helpful counseling point.")

# Option III: Highlights that risks are not understood in adolescents and they should not vape.
# This is the core educational message, differentiating the son's situation from the mother's.
option_III_is_correct = True
print("Statement III provides the essential rationale for why the adolescent should stop vaping.")

# The best answer combines the rationale for quitting (III) with a first-line treatment plan (II).
# The final answer is the combination of these two statements.
final_combination = ["II", "III"]

print(f"The best combination of counseling points is {', '.join(final_combination)}.")

# Restore stdout
sys.stdout = old_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the thinking process
print(output)
# The final answer choice is J, which corresponds to the combination of II and III.
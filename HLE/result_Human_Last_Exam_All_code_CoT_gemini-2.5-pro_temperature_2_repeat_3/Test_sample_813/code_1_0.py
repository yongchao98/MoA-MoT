import sys
import io

# Store the original stdout
original_stdout = sys.stdout
# Create a new StringIO object
captured_output = io.StringIO()
# Redirect stdout
sys.stdout = captured_output

# --- Start of the analysis code ---
age = 59
work_history_years = 30

# Step 1: Analyze the initial symptoms and patient history.
# The combination of symptoms (agitation, hypersexuality, etc.) and a family history of mood disorders
# points to a diagnosis of a manic episode, likely related to Bipolar Disorder.
analysis_step_1 = "Patient's initial symptoms are characteristic of a manic episode."

# Step 2: Identify the likely intervention.
# Lithium is a first-line medication used to treat Bipolar Disorder.
analysis_step_2 = "A common and likely prescription for these symptoms is Lithium."

# Step 3: Connect the intervention to a known side effect.
# A major side effect of Lithium is drug-induced hypothyroidism.
analysis_step_3 = "A well-documented side effect of Lithium is hypothyroidism."

# Step 4: Connect the side effect to the final symptom.
# A classic symptom of hypothyroidism is decreased libido or sexual dysfunction.
analysis_step_4 = "Hypothyroidism is a known cause of decreased sexual interest."

# Step 5: Conclude the causal chain.
conclusion = "The most logical root cause is that Lithium, prescribed for the manic episode, induced hypothyroidism, which in turn caused the sexual dysfunction."
final_answer_choice = "A"

print("Step-by-Step Reasoning:")
print(f"Patient Age: {age}")
print(f"Work History: {work_history_years} years in metal smelting")
print(f"1. {analysis_step_1}")
print(f"2. {analysis_step_2}")
print(f"3. {analysis_step_3}")
print(f"4. {analysis_step_4}")
print("\nConclusion:")
print(f"The causal chain strongly supports that the underlying root cause is 'Lithium induced hypothyroidism'.")
# --- End of the analysis code ---

# Restore the original stdout
sys.stdout = original_stdout
# Get the content of the captured output
output = captured_output.getvalue()

# Print the captured output to the actual console
print(output, end="")

# Final answer in the required format
print(f"<<<{final_answer_choice}>>>")
import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a new string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

def solve_opioid_tapering_case():
    """
    This function analyzes the clinical scenario and prints a step-by-step solution.
    """

    # 1. Define the detailed step-by-step reasoning.
    explanation = """Here is the step-by-step analysis to determine the best course of action for the patient:

1.  **Analyze the Patient's Situation:** The patient has a complex medical history (remission from lymphoma) and is on high-dose opioids with documented challenges in tapering. This situation strongly suggests significant physical dependence and a high likelihood of an underlying Opioid Use Disorder (OUD), which requires a comprehensive and specialized approach.

2.  **Evaluate the Options:**
    *   **Statement I (Simple Taper):** Incorrect. A simple taper of the current opioid is unlikely to succeed, as the patient is already struggling. It does not address the potential underlying OUD.
    *   **Statement II (Methadone):** Correct. Methadone is a long-acting opioid and a first-line, evidence-based treatment for both OUD and complex chronic pain. It is a very appropriate option for the multidisciplinary team to consider.
    *   **Statement III (Rapid Taper):** Incorrect. A rapid taper from high-dose opioids is dangerous, can cause severe withdrawal, and increases the risk of relapse and overdose. This is not a safe or effective strategy.
    *   **Statement IV (Multidisciplinary Consultation):** Correct. A consultation involving pain management and psychiatry/addiction medicine is the absolute standard of care for a complex case like this. It is essential for a safe and effective plan.
    *   **Statement V (Buprenorphine-Naloxone):** Correct. Buprenorphine-naloxone is a first-line treatment for OUD with a superior safety profile to full agonists. It effectively manages withdrawal and cravings and directly addresses the patient's specific question.

3.  **Synthesize the Best Plan:** The ideal management plan involves a multidisciplinary team (IV) to conduct a thorough assessment. This team would then consider the premier, evidence-based Medication-Assisted Treatment (MAT) options, which are methadone (II) and buprenorphine-naloxone (V). The selection of all three statements represents the most comprehensive and appropriate approach.

4.  **Conclusion:** The best combination of statements describing the optimal approach is II, IV, and V."""

    # 2. Identify the selected statements for the "equation" and the final answer choice.
    selected_statements_numerals = ["II", "IV", "V"]
    final_answer_choice = "H"

    # 3. Print the full explanation.
    print(explanation)
    
    # 4. Fulfill the request to output the "equation" of the selected numbers.
    final_equation = " + ".join(selected_statements_numerals)
    print("\nThe best plan incorporates the principles from statements:")
    print(final_equation)

    # 5. Print the final answer in the required format.
    print(f"\n<<<{final_answer_choice}>>>")

# Execute the function to generate the output
solve_opioid_tapering_case()

# Restore original stdout
sys.stdout = original_stdout
# Get the content from the buffer
final_output = output_buffer.getvalue()

# Print the final output to the actual console
print(final_output)
import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_clinical_scenario():
    """
    Analyzes the clinical scenario and prints the rationale and final answer.
    """
    # Define the chosen options based on clinical reasoning
    prioritized_options = {
        "I": "Counsel patient on stopping cannabis.",
        "III": "Order a urine drug test.",
        "IV": "Prescribe melatonin for insomnia."
    }
    
    final_answer_letter = "L"

    # Print the detailed rationale
    print("Rationale for the prioritized treatment options:")
    print("-" * 40)
    
    print("Option I: Counsel patient on stopping cannabis.")
    print("This is a critical first step. The patient's heavy cannabis use is a major confounding variable that likely worsens his anxiety and insomnia, and interferes with the effectiveness of his medications.")
    print("\n")
    
    print("Option III: Order a urine drug test.")
    print("This is essential for gathering objective data. It confirms the self-reported cannabis use and screens for other potential substances, which is standard procedure given the patient's history of multiple substance use disorders.")
    print("\n")

    print("Option IV: Prescribe melatonin for insomnia.")
    print("This directly addresses a severe symptom (insomnia) with a safe, low-risk, non-addictive intervention. Providing some symptomatic relief can help build the therapeutic alliance needed to tackle the more complex underlying issues.")
    print("-" * 40)
    
    # Print the selected option numbers as requested
    print("The three prioritized options are I, III, and IV.")
    
    # Print the final answer in the specified format
    print(f"<<<{final_answer_letter}>>>")

# Execute the function to generate the output
solve_clinical_scenario()

# Restore stdout and print the captured output to the console
sys.stdout = old_stdout
print(captured_output.getvalue())
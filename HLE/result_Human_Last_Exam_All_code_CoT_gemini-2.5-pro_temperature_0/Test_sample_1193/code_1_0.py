import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring options based on key patient data.
    """
    patient_data = {
        "procedure": "Whipple procedure",
        "time_post_op_days": 29,
        "oxygen_saturation_percent": 82,
        "supplemental_oxygen_L": 3,
        "respiratory_effort": "gasping for air",
        "lung_sounds": "bilateral crackles"
    }

    diagnoses = {
        "A": {"name": "Acute blood transfusion reaction", "reasoning": "This typically occurs within hours to a day of a transfusion. The patient's timeline of 29 days makes this highly unlikely."},
        "B": {"name": "Iodine-related reaction", "reasoning": "This is a reaction to contrast dye. There is no mention of a recent contrast study, and the presentation does not fit."},
        "C": {"name": "Sensitivity reaction", "reasoning": "This term is too vague and does not explain the specific, severe syndrome of bilateral crackles and severe hypoxemia."},
        "D": {"name": "Sepsis", "reasoning": "The patient underwent major surgery, a risk factor for infection. The timeline of 29 days is sufficient for a post-operative infection to develop. Sepsis is the most common cause of Acute Respiratory Distress Syndrome (ARDS), which perfectly matches the clinical picture: severe hypoxemia (O2 sat: 82%), respiratory distress, and bilateral crackles (pulmonary edema)."},
        "E": {"name": "Myocyte necrosis", "reasoning": "This does not directly cause this respiratory presentation. While a severe heart attack could cause pulmonary edema, sepsis-induced ARDS is a more common complication in this post-operative setting."},
        "F": {"name": "Respiratory deconditioning", "reasoning": "This is a gradual weakening and would not cause such acute, severe hypoxemia or the bilateral crackles indicating fluid in the lungs."},
        "G": {"name": "Lung exhaustion", "reasoning": "This is not a standard medical diagnosis and does not explain the underlying cause of the fluid in the lungs."},
        "H": {"name": "Air pollution sensitivity", "reasoning": "This is highly unlikely to be the cause of such a severe clinical deterioration in a post-operative patient."}
    }

    print("Analyzing the patient's case based on the following data:")
    print(f"- Time since surgery: {patient_data['time_post_op_days']} days")
    print(f"- Oxygen Level: {patient_data['oxygen_saturation_percent']}% on {patient_data['supplemental_oxygen_L']}L of oxygen")
    print(f"- Physical Exam: '{patient_data['respiratory_effort']}' and '{patient_data['lung_sounds']}'")
    print("-" * 30)
    print("Evaluation of Potential Causes:")
    
    # The most likely diagnosis based on clinical reasoning
    most_likely_choice = "D"
    
    for choice, details in diagnoses.items():
        print(f"Choice {choice} ({details['name']}): {details['reasoning']}")
        if choice == most_likely_choice:
            print(f"--> Conclusion: This is the most plausible diagnosis.\n")
        else:
            print(f"--> Conclusion: This is less likely.\n")

    print("-" * 30)
    print(f"Final Answer: The most likely cause of the patient's hypoxemia is Sepsis (leading to ARDS).")


solve_medical_case()

# Restore stdout and print the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()
print(output)
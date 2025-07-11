import sys
import io

# Redirect stdout to capture the output of the function
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def diagnose_patient():
    """
    Analyzes clinical findings to determine the most likely disease by scoring diagnoses.
    """
    # Key clinical findings extracted from the case vignette
    findings = {
        "Immunosuppression (steroids)": 10,
        "Pulmonary Nodules": 8,
        "Acute Pneumonia & Sepsis": 9,
        "Cutaneous Lesions": 8,
        "CNS Symptoms (confusion)": 7,
        "Failure of Aminoglycoside Therapy": 10,
        "Underlying Chronic Polyarthritis": 6,
        "Relevant History (smoking, asbestos risk)": 5
    }

    # Differential diagnoses and how well they match the findings (on a scale of 0-10)
    diagnoses = {
        "Nocardiosis": {
            "Immunosuppression (steroids)": 10,  # Classic host factor
            "Pulmonary Nodules": 9,             # Common presentation
            "Acute Pneumonia & Sepsis": 9,      # Common presentation
            "Cutaneous Lesions": 9,             # Common with dissemination
            "CNS Symptoms (confusion)": 8,      # Common with dissemination (brain abscess)
            "Failure of Aminoglycoside Therapy": 10, # Nocardia is intrinsically resistant
            "Underlying Chronic Polyarthritis": 2, # Can cause septic arthritis, but not chronic inflammatory arthritis
            "Relevant History (smoking, asbestos risk)": 4 # Chronic lung damage is a risk factor
        },
        "Granulomatosis with Polyangiitis (GPA)": {
            "Immunosuppression (steroids)": 10, # GPA is the reason for steroid use
            "Pulmonary Nodules": 10,            # Classic feature
            "Acute Pneumonia & Sepsis": 1,      # GPA is inflammation, not infection
            "Cutaneous Lesions": 7,             # Vasculitic lesions are possible
            "CNS Symptoms (confusion)": 8,      # CNS vasculitis is possible
            "Failure of Aminoglycoside Therapy": 0, # Not an infection, so not applicable
            "Underlying Chronic Polyarthritis": 9, # Classic feature
            "Relevant History (smoking, asbestos risk)": 7 # Known associations
        },
        "Tuberculosis (TB)": {
            "Immunosuppression (steroids)": 8,
            "Pulmonary Nodules": 8,
            "Acute Pneumonia & Sepsis": 7,
            "Cutaneous Lesions": 5,
            "CNS Symptoms (confusion)": 6,
            "Failure of Aminoglycoside Therapy": 2, # Aminoglycosides (e.g., streptomycin) ARE used to treat TB
            "Underlying Chronic Polyarthritis": 3,
            "Relevant History (smoking, asbestos risk)": 5
        }
    }

    scores = {disease: 0 for disease in diagnoses}

    print("Diagnostic Reasoning:")
    print("="*25)
    print("Scoring potential diagnoses based on weighted clinical findings.")
    print("The final diagnosis should explain the entire clinical picture, especially the terminal infectious illness.\n")

    # This loop simulates creating a 'final equation' by showing each component's contribution
    for finding, weight in findings.items():
        print(f"Factor: '{finding}' | Weight: {weight}")
        equation_parts = []
        for disease in diagnoses:
            match_score = diagnoses[disease].get(finding, 0)
            calculated_score = match_score * weight
            scores[disease] += calculated_score
            equation_parts.append(f"{disease}[+{calculated_score}]")
        print("  " + ", ".join(equation_parts))
        print("-" * 20)

    # Note on GPA: GPA is the likely underlying condition that led to immunosuppression, making Nocardiosis possible.
    # The question "What disease did this person likely experience?" points towards the entire process,
    # culminating in the specific, identifiable cause of death.
    final_diagnosis = "Nocardiosis (superimposed on an underlying disease like GPA)"

    print("\nFinal Score Calculation:")
    print("="*25)
    print("Score [GPA]: This autoimmune disease explains the chronic initial symptoms (arthritis, nodules) and the need for steroids, but not the acute infectious death. Score =", scores['Granulomatosis with Polyangiitis (GPA)'])
    print("Score [TB]: Fails to explain the resistance to Aminoglycoside therapy. Score =", scores['Tuberculosis (TB)'])
    print("Score [Nocardiosis]: Best explains the opportunistic infection in an immunocompromised host, including pulmonary, cutaneous, and CNS signs, plus the specific antibiotic resistance. Score =", scores['Nocardiosis'])
    print("\nConclusion:")
    print(f"The most likely diagnosis for the terminal septic shock is {final_diagnosis.split('(')[0].strip()}.")


diagnose_patient()

# Restore stdout and print the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()
print(output)
<<<Nocardiosis>>>
import sys
import io

# Capture the original stdout
original_stdout = sys.stdout
# Create a new StringIO object
new_stdout = io.StringIO()
# Redirect stdout
sys.stdout = new_stdout

def diagnose_patient_case():
    """
    Analyzes the clinical case and prints the most likely diagnosis.
    """
    # Key factors from the clinical case leading to the diagnosis
    patient_profile = "62-year-old man, smoker, history of shipbuilding (asbestos exposure)."
    initial_symptoms = "Systemic illness: polyarthritis, fatigue, confusion, bruising, and multiple pulmonary nodules."
    treatment_course = "Initial improvement with steroids, indicating an inflammatory/autoimmune process."
    final_event = "Development of a severe, fatal opportunistic infection resistant to aminoglycosides after steroid treatment and travel."

    # Synthesis of the clinical data
    diagnosis = "Granulomatosis with Polyangiitis (GPA), also known as Wegener's Granulomatosis."
    reasoning = (
        "This systemic autoimmune vasculitis best explains the full constellation of symptoms. "
        "It accounts for the pulmonary nodules and the multi-organ inflammatory signs (joints, skin, neurological). "
        "Crucially, the disease and its treatment with steroids lead to severe immunosuppression, "
        "making the patient vulnerable to the fatal opportunistic infection he contracted."
    )

    print(f"The patient's clinical presentation, including pulmonary nodules, systemic vasculitis, "
          f"and fatal opportunistic infection following steroid therapy, is most characteristic "
          f"of the following disease:")
    print(f"\nLikely Diagnosis: {diagnosis}")
    print(f"\nReasoning: {reasoning}")

diagnose_patient_case()

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output = new_stdout.getvalue()

# Print the captured output to the actual console
print(output)

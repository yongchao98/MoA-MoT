import sys
import io

# A simple function to print the analysis to the console.
def solve_clinical_case():
    """
    This script analyzes a clinical vignette to determine the correct diagnostic maneuver.
    It breaks down the problem into logical steps:
    1. Parsing the case details.
    2. Identifying the probable diagnosis.
    3. Evaluating physical exam maneuvers based on anatomy and function.
    4. Concluding with the correct answer.
    """
    
    # Step 1: Define the key clinical findings from the vignette.
    symptoms = "Pain in the lower right extremity, L4-S1 distribution."
    imaging = "X-ray imaging is unremarkable."
    exam_position = "Patient in left decubitus position, testing the extended right leg."
    probable_diagnosis = "Piriformis Syndrome"

    print("--- Clinical Case Analysis ---")
    print(f"Symptoms: {symptoms}")
    print(f"Imaging Findings: {imaging}")
    print(f"Physical Exam Setup: {exam_position}")
    print("-" * 30)

    # Step 2: State the working diagnosis based on the findings.
    print(f"Probable Diagnosis: {probable_diagnosis}")
    print("Reasoning: The symptoms are consistent with sciatica. Unremarkable X-rays suggest a soft tissue cause, such as nerve entrapment by the piriformis muscle.")
    print("-" * 30)

    # Step 3: Analyze the relevant anatomy and choose the confirmatory test.
    print("Analyzing the required action:")
    print("The piriformis muscle is a primary EXTERNAL rotator of the extended hip.")
    print("To confirm piriformis syndrome, one must provoke the pain by stressing the muscle.")
    print("This can be done by:")
    print("  a) Stretching the muscle -> Forcing the hip into INTERNAL rotation.")
    print("  b) Contracting the muscle -> Resisting the patient's EXTERNAL rotation.")
    print("\nBoth are valid tests. The maneuver of passively stretching the piriformis is a classic diagnostic sign (Freiberg's sign).")
    
    # Step 4: Map the correct action to the answer choices.
    answer_choices = {
        'A': 'Abduction',
        'B': 'Adduction',
        'C': 'Internal Rotation',
        'D': 'External Rotation',
        'E': 'Flexion',
        'F': 'Extension'
    }
    
    correct_action = 'Internal Rotation'
    correct_choice_letter = 'C'

    print(f"\nThe action that stretches the piriformis and confirms the diagnosis is '{correct_action}'.")
    print(f"This corresponds to answer choice: {correct_choice_letter}")
    print("-" * 30)


# Run the solver
solve_clinical_case()
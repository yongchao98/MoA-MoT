import sys

# This script simulates the clinical reasoning for managing an agitated patient.
# WARNING: This is not medical advice. Consult a qualified healthcare professional for medical emergencies.

def evaluate_agitation_options():
    """
    Analyzes the best next step for an agitated patient who has failed initial treatment.
    """
    # Patient Information
    patient_status = "Acutely agitated and violent, punched a physician."
    unknown_history = True
    initial_treatment = "5mg of Zyprexa (olanzapine) IM"
    initial_outcome = "No improvement"

    print("--- Patient Case Analysis ---")
    print(f"Status: {patient_status}")
    print(f"History: Unknown medical and allergy history: {unknown_history}")
    print(f"Initial Treatment: {initial_treatment}")
    print(f"Outcome: {initial_outcome}\n")

    print("--- Evaluating Next Step Options ---")

    # Option A
    print("\n[Option A] 2mg IV lorazepam")
    print("Analysis: Requires IV access, which is difficult and potentially dangerous to obtain in a combative patient. IM route is generally safer for staff in this scenario.")

    # Option B
    print("\n[Option B] 2mg IM lorazepam + 5mg olanzapine IM")
    print("Analysis: This adds a benzodiazepine (lorazepam) to the antipsychotic (olanzapine), a common strategy for a synergistic effect. The total olanzapine dose becomes 10mg (5mg initial + 5mg new), which is a standard dose. This is a balanced approach to escalation.")

    # Option C
    print("\n[Option C] Verbal de-escalation before any pharmacologic intervention")
    print("Analysis: The patient has already become physically violent. The window for verbal de-escalation as the sole intervention has passed. Safety is the priority.")

    # Option D
    print("\n[Option D] 10mg IM olanzapine")
    print("Analysis: This escalates a single agent that has already failed at a lower dose. The total olanzapine dose would be 15mg (5mg initial + 10mg new), a high dose for a patient with an unknown history.")

    # Option E
    print("\n[Option E] 10mg IM olanzapine + 2mg IM lorazepam")
    print("Analysis: This is a highly sedating combination. The total olanzapine dose of 15mg plus lorazepam carries a significant risk of over-sedation and cardiorespiratory depression, which is dangerous in a patient with an unknown background.")

    print("\n--- Conclusion ---")
    print("Option B is the most appropriate next step. It provides effective combination therapy by adding a second agent after the first was ineffective, using standard doses that balance efficacy with safety.")
    
    print("\nFinal Recommended Intervention:")
    # The user prompt requested printing the numbers in the final equation.
    # The chosen intervention is B: 2mg IM lorazepam + 5mg olanzapine IM
    lorazepam_dose = 2
    olanzapine_dose = 5
    print(f"{lorazepam_dose}mg IM lorazepam + {olanzapine_dose}mg olanzapine IM")

# Execute the analysis
evaluate_agitation_options()
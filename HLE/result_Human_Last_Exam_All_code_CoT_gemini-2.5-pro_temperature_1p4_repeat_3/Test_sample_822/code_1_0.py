import sys
import io

# Redirect stdout to capture the output for the final answer tag
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def provide_diagnosis():
    """
    Analyzes the provided clinical case and prints the step-by-step reasoning
    and the final diagnosis.
    """
    print("Analyzing the clinical case to determine the most likely disease...")
    print("-" * 20)

    print("Key Clues from the Case:")
    print("1. Patient Profile: 62-year-old male, smoker, worked in shipbuilding (asbestos exposure).")
    print("2. Initial Symptoms: Systemic inflammation with polyarthritis (wrists, ankles, elbows), fatigue, confusion, bruising.")
    print("3. Key Finding: Chest X-ray showing multiple pulmonary nodules.")
    print("4. Clinical Course: Started on steroids (immunosuppression), then developed a severe opportunistic infection (fever, cough, cutaneous lesions) after travel.")
    print("5. Treatment Failure: The infection did not respond to Aminoglycoside therapy.")
    print("6. Outcome: Death from septic shock, indicating an overwhelming infection in a compromised host.")
    print("-" * 20)

    print("Diagnostic Reasoning:")
    print("The constellation of symptoms points strongly to a systemic autoimmune vasculitis. Let's evaluate the best fit.")
    print("- Wegener's Granulomatosis, now known as Granulomatosis with Polyangiitis (GPA), is a vasculitis that classically presents with a triad affecting the lungs, kidneys, and upper respiratory tract.")
    print("- The patient's presentation of multiple pulmonary nodules and polyarthritis is highly characteristic of GPA.")
    print("- Treatment for GPA involves steroids and other immunosuppressants, which explains why the patient was started on this therapy.")
    print("- A major complication of GPA, especially when treated with immunosuppressants, is a severe opportunistic infection. The patient's final illness (pulmonary and cutaneous lesions, resistance to aminoglycosides) is a classic picture of disseminated Nocardiosis, a well-known and often fatal complication in patients with GPA.")
    print("-" * 20)

    print("Conclusion:")
    final_diagnosis = "Granulomatosis with Polyangiitis (GPA)"
    print(f"The patient's entire clinical course, from the initial presentation of pulmonary nodules and arthritis to the fatal opportunistic infection following immunosuppressive therapy, makes {final_diagnosis} the most likely diagnosis.")

provide_diagnosis()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the captured output to the console
print(output)

# Extract the final answer for the tag
final_answer_line = [line for line in output.split('\n') if "most likely diagnosis" in line]
# A bit of parsing to be robust
final_answer = final_answer_line[0].split("makes ")[1].split(" the most")[0] if final_answer_line else "Granulomatosis with Polyangiitis (GPA)"
print(f"<<<{final_answer}>>>")
import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer to format the final output correctly
sys.stdout = captured_output = io.StringIO()

def analyze_dermatology_case():
    """
    This function analyzes the provided clinical case and determines the best next diagnostic step.
    It prints the reasoning based on key facts from the case.
    """
    # Key details from the case
    patient_age = 43
    rash_location = "posterior border of both axillary folds"
    spared_area = "axillary vaults"
    key_history_fact = "started a weight-loss workout program 3 weeks ago"

    print("Analyzing the clinical case to determine the best next diagnostic step:")
    print("="*60)

    # Step 1: Correlate patient history with symptoms
    print(f"1. History Analysis:")
    print(f"   - A {patient_age}-year-old man presents with an intensely itchy rash.")
    print(f"   - The patient {key_history_fact}, which likely involves new workout clothes and increased perspiration.")
    print(f"   - This timing creates a strong suspicion for an allergic reaction to his new clothing (textile dermatitis).")
    print("")

    # Step 2: Analyze the physical exam findings
    print("2. Physical Exam Analysis:")
    print(f"   - The rash is located on the '{rash_location}'.")
    print(f"   - Crucially, the '{spared_area}' are spared.")
    print(f"   - This specific distribution points to an external agent rubbing against the skin folds (like a t-shirt) rather than a product applied directly into the armpit (like deodorant).")
    print("")

    # Step 3: Evaluate the diagnostic options provided
    print("3. Evaluating the Answer Choices:")
    print("   - A. Skin biopsy: Too invasive for a primary diagnostic step in this classic presentation.")
    print("   - B. KOH preparation: Tests for fungus. Less likely given the specific distribution and history.")
    print("   - C. Topical steroid: This is a form of treatment, not a diagnostic test.")
    print("   - D. Patch test: This is the definitive test (gold standard) to identify the specific allergen causing Allergic Contact Dermatitis (ACD).")
    print("")

    # Step 4: Formulate a conclusion
    print("4. Conclusion:")
    print("   The clinical picture strongly suggests Allergic Contact Dermatitis from textiles.")
    print("   The best next step to confirm this diagnosis and identify the causative agent is a patch test.")
    print("="*60)
    print("The correct choice is D.")

# Run the analysis
analyze_dermatology_case()

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the formatted output as requested
print(output)
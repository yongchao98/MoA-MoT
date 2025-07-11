import sys
import io

# Redirect stdout to a string buffer to capture all prints
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_case():
    """
    Analyzes the patient's clinical information to determine the most likely diagnosis.
    """
    # --- Patient Data & Risk Factor Calculation ---
    patient_age = 53
    pain_level = 10
    cigarette_packs_per_day = 2
    smoking_duration_years = 20

    # The problem requires showing an equation with numbers from the prompt.
    # Calculating pack-years is a clinically relevant equation.
    pack_years = cigarette_packs_per_day * smoking_duration_years

    print("Step-by-Step Diagnostic Analysis")
    print("=" * 35)

    # --- Step 1: Analyze Patient Risk Factors ---
    print("\n1. Assessing Key Patient Risk Factors:")
    print(f" - Patient is {patient_age} years old.")
    print(f" - Significant history of alcohol use disorder.")
    print(f" - Major smoking history is quantified using the following equation:")
    print(f"   {cigarette_packs_per_day} packs/day * {smoking_duration_years} years = {pack_years} pack-years")
    print("   (This represents a very high risk for squamous cell carcinoma).")

    # --- Step 2: Correlate Symptoms and Findings ---
    print("\n2. Correlating Symptoms with Diagnostic Results:")
    print(f" - Symptoms: Severe chest pain ({pain_level}/10) and odynophagia (painful swallowing).")
    print(f" - Imaging: Shows esophageal wall thickening and narrowing, suggesting an infiltrative process.")
    print(f" - Endoscopy: Negative for mucosal lesions (ulcers, plaques). This is a critical finding.")

    # --- Step 3: Evaluating the Differential Diagnosis ---
    print("\n3. Evaluating Potential Diagnoses:")
    print(" - A/E (Infectious Esophagitis): Unlikely, as endoscopy would typically show ulcers or plaques.")
    print(" - D (GERD): Unlikely, given the negative endoscopy and imaging findings of wall thickening.")
    print(" - B (Esophageal Adenocarcinoma): Less likely. The patient's risk profile (smoking, alcohol) strongly favors SCC.")
    print(" - C (Esophageal Squamous Cell Carcinoma): Most likely. The powerful risk factors ({pack_years} pack-years smoking, alcohol) combined with symptoms and imaging are highly suggestive. The negative endoscopy can be explained by a submucosal or infiltrating tumor that has thickened the wall without yet ulcerating the surface.")

    # --- Step 4: Conclusion ---
    print("\nConclusion:")
    print("The combination of overwhelming risk factors and findings consistent with an infiltrating tumor makes Esophageal Squamous Cell Carcinoma the most probable diagnosis.")

# Execute the analysis
analyze_case()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())

# Final Answer in the required format
print("<<<C>>>")
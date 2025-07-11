import sys
import io

# Redirect stdout to capture prints
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_medical_case():
    """
    Analyzes a clinical case to determine the most likely diagnosis based on a scoring system.
    """
    # --- Patient Data from the Case ---
    age = 53
    pain_intensity = 10
    smoke_packs_per_day = 2
    smoke_years = 20
    pack_years = smoke_packs_per_day * smoke_years

    # --- Initialize scores for each diagnosis ---
    # A. Streptococcal esophagitis
    # B. Esophageal adenocarcinoma
    # C. Esophageal squamous cell carcinoma
    # D. GERD
    # E. Herpes esophagitis
    diagnoses = {
        'A': 0, 'B': 0, 'C': 0, 'D': 0, 'E': 0
    }
    
    # --- Scoring based on clinical evidence ---

    # 1. Major Risk Factors (Smoking & Alcohol)
    # The patient's 40 pack-year smoking history and alcohol use are massive risk factors for SCC.
    diagnoses['C'] += 10

    # 2. Imaging Findings (Wall Thickening & Narrowing)
    # This points to an infiltrative process, common in cancer.
    diagnoses['B'] += 5
    diagnoses['C'] += 5

    # 3. Endoscopy Result (Normal Mucosa)
    # This is the key finding. It rules out conditions that cause visible surface changes.
    # It strongly supports an infiltrative tumor that grows beneath the surface.
    diagnoses['A'] -= 10 # Would have plaques
    diagnoses['D'] -= 10 # Would have redness/erosions
    diagnoses['E'] -= 10 # Would have ulcers
    diagnoses['B'] -= 5  # Usually presents as a visible mass
    diagnoses['C'] += 10 # Classic for infiltrative SCC

    # --- Output the Analysis ---
    
    print("### Clinical Case Analysis ###")
    print("\nStep 1: Evaluating key patient data points for the equation.")
    print(f"Patient's age: {age}")
    print(f"Pain intensity: {pain_intensity}")
    print(f"Smoking history: {smoke_packs_per_day} packs/day for {smoke_years} years = {pack_years} pack-years")

    print("\nStep 2: Calculating diagnostic scores based on the evidence.")
    print("The final 'equation' is the sum of scores for each diagnosis:")
    
    final_equation_output = []
    for diagnosis, score in diagnoses.items():
        final_equation_output.append(f"Diagnosis {diagnosis} Final Score = {score}")
    print("\n".join(final_equation_output))
    
    most_likely = max(diagnoses, key=diagnoses.get)
    
    print("\nStep 3: Conclusion")
    print(f"\nThe most likely diagnosis is C (Esophageal squamous cell carcinoma).")
    print("Reasoning: The patient's major risk factors (heavy smoking, alcohol use) combined with imaging showing wall thickening and a NORMAL endoscopy is a classic presentation for an infiltrative squamous cell carcinoma, which grows within the esophageal wall without breaking through the surface lining.")

solve_medical_case()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
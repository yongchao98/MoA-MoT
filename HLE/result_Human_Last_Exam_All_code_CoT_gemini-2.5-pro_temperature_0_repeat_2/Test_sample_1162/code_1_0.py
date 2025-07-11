import sys
import io

# Redirect stdout to capture the output for the final answer format
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_diagnosis():
    """
    Analyzes clinical findings to determine the most likely diagnosis using a weighted scoring system.
    """
    # 1. Define the patient's key clinical findings from the vignette
    patient_findings = {
        "aniridia": True,
        "hypertension": True,
        "pelvic_mass": True,
        "developmental_delay": True,
        "conjunctival_pallor": True, # General sign of anemia
        "lethargy": True # General sign
    }

    # 2. Define diagnostic knowledge base with weighted associations.
    # A high weight (e.g., 10) indicates a very strong or classic association.
    # A pelvic mass is treated as a possible abdominal/renal mass for scoring.
    diagnostic_associations = {
        "A. Germ cell tumor": {"pelvic_mass": 5},
        "B. Astrocytoma": {"neurological_deficits": 10},
        "C. Neuroblastoma": {"abdominal_mass": 4, "hypertension": 3, "opsoclonus_myoclonus": 10},
        "D. Nephroblastoma": {"abdominal_mass": 5, "hypertension": 3, "aniridia": 10, "developmental_delay": 2},
        "E. Ewing sarcoma": {"bone_or_soft_tissue_mass": 5}
    }

    # Map the patient's pelvic mass to a more general category for scoring
    if patient_findings.get("pelvic_mass"):
        patient_findings["abdominal_mass"] = True
        patient_findings["bone_or_soft_tissue_mass"] = True

    # 3. Calculate scores and find the best match
    best_diagnosis = ""
    max_score = -1

    for diagnosis, associations in diagnostic_associations.items():
        current_score = 0
        for finding, present in patient_findings.items():
            if present and finding in associations:
                current_score += associations[finding]
        
        if current_score > max_score:
            max_score = current_score
            best_diagnosis = diagnosis

    # 4. Output the reasoning and the final conclusion
    print("Analyzing the clinical vignette to find the most probable diagnosis...")
    print("-" * 60)
    print(f"The diagnosis with the highest match score is: {best_diagnosis}")
    print("\nThis conclusion is based on the following key findings:")
    
    # Build and print the final equation for the best diagnosis
    equation_parts = []
    final_score = 0
    for finding, score in diagnostic_associations[best_diagnosis].items():
        # Check if the finding is present in the patient's case
        if patient_findings.get(finding):
            # Add the finding and its score to the equation string
            equation_parts.append(f"{finding.replace('_', ' ').title()}[{score}]")
            final_score += score
    
    # Print the final equation showing each number
    final_equation_str = " + ".join(equation_parts) + f" = {final_score}"
    print("\nDiagnostic Score Equation:")
    print(final_equation_str)

    print("\nExplanation:")
    print("The patient's symptoms strongly suggest WAGR syndrome, which includes:")
    print(" - (W)ilms Tumor: A kidney tumor, also known as Nephroblastoma. This presents as an abdominal/pelvic mass and can cause hypertension.")
    print(" - (A)niridia: Absence of the iris. This is a very specific finding and is present in the patient.")
    print(" - (G)enitourinary anomalies.")
    print(" - (R)etardation: The patient's delayed speech fits this criterion.")
    print("\nTherefore, the most likely diagnosis is Nephroblastoma (Wilms tumor).")

solve_diagnosis()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()
print(output)

# Final answer in the required format
print("<<<D>>>")
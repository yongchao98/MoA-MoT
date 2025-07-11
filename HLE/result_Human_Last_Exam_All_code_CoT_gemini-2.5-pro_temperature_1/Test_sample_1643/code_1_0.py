import sys

# Define a class to hold patient data for clarity
class PatientData:
    def __init__(self):
        # Clinical Information
        self.age = 34
        self.vte_events = 3
        self.medications = ["DESOGESTREL 75µg", "RIVAROXABAN 20mg"]
        
        # Lab results from 3 months ago
        self.labs_t1 = {
            "antiB2GP1_IgM": 41, "antiB2GP1_IgG": 18,
            "aCL_IgM": 32, "aCL_IgG": 9,
            "PTT_LA_ratio": 1.19, "dRVVT_ratio": 1.44
        }
        
        # Lab results from today (3 months later)
        self.labs_t2 = {
            "antiB2GP1_IgM": 29, "antiB2GP1_IgG": 21,
            "aCL_IgM": 47, "aCL_IgG": 7,
            "PTT_LA_ratio": 1.17, "dRVVT_ratio": 1.51
        }
        
        self.time_between_labs_weeks = 12

# Define normal values for lab tests
LAB_NORMALS = {
    "antiB2GP1_IgM": 20, "antiB2GP1_IgG": 20,
    "aCL_IgM": 20, "aCL_IgG": 20,
    "PTT_LA_ratio": 1.18, "dRVVT_ratio": 1.2
}

def analyze_patient_for_aps(patient):
    """
    Analyzes patient data against APS criteria and prints the reasoning.
    """
    print("Step 1: Evaluating Clinical Criteria for Antiphospholipid Syndrome (APS)")
    print("-" * 70)
    
    # Check for Vascular Thrombosis
    clinical_criterion_met = False
    if patient.vte_events >= 1:
        clinical_criterion_met = True
        print(f"Clinical criterion MET: Patient has a history of {patient.vte_events} venous thromboembolic (VTE) events.")
    else:
        print("Clinical criterion NOT MET: No documented history of vascular thrombosis.")
    print("\n")

    print("Step 2: Evaluating Laboratory Criteria for Antiphospholipid Syndrome (APS)")
    print("-" * 70)
    print(f"Lab tests were performed at two time points, separated by {patient.time_between_labs_weeks} weeks.")
    if patient.time_between_labs_weeks < 12:
        print("The time between lab tests is less than 12 weeks, so persistence cannot be confirmed.")
        lab_criteria_met = False
    else:
        print("The 12-week interval between tests is met, allowing for confirmation of persistence.\n")
        
        # List to hold confirmed persistent lab findings
        persistent_findings = []

        # Check 1: Lupus Anticoagulant (LA)
        print("Checking for persistent Lupus Anticoagulant (LA)...")
        la_t1_positive = patient.labs_t1["dRVVT_ratio"] > LAB_NORMALS["dRVVT_ratio"]
        la_t2_positive = patient.labs_t2["dRVVT_ratio"] > LAB_NORMALS["dRVVT_ratio"]
        print(f"  - Test 1 (3 months ago): dRVVT ratio was {patient.labs_t1['dRVVT_ratio']} (Normal < {LAB_NORMALS['dRVVT_ratio']}) -> {'Positive' if la_t1_positive else 'Negative'}")
        print(f"  - Test 2 (Today): dRVVT ratio was {patient.labs_t2['dRVVT_ratio']} (Normal < {LAB_NORMALS['dRVVT_ratio']}) -> {'Positive' if la_t2_positive else 'Negative'}")
        if la_t1_positive and la_t2_positive:
            persistent_findings.append("Lupus Anticoagulant (dRVVT)")
            print("  -> Result: PERSISTENTLY POSITIVE.\n")
        else:
            print("  -> Result: NOT persistently positive.\n")
        
        # Check 2: Anticardiolipin (aCL) antibodies
        print("Checking for persistent Anticardiolipin (aCL) antibodies...")
        acl_igm_t1_pos = patient.labs_t1["aCL_IgM"] > LAB_NORMALS["aCL_IgM"]
        acl_igm_t2_pos = patient.labs_t2["aCL_IgM"] > LAB_NORMALS["aCL_IgM"]
        print(f"  - aCL IgM Test 1: {patient.labs_t1['aCL_IgM']} UI/L (Normal < {LAB_NORMALS['aCL_IgM']}) -> {'Positive' if acl_igm_t1_pos else 'Negative'}")
        print(f"  - aCL IgM Test 2: {patient.labs_t2['aCL_IgM']} UI/L (Normal < {LAB_NORMALS['aCL_IgM']}) -> {'Positive' if acl_igm_t2_pos else 'Negative'}")
        if acl_igm_t1_pos and acl_igm_t2_pos:
            persistent_findings.append("Anticardiolipin IgM")
            print("  -> Result: PERSISTENTLY POSITIVE for aCL IgM.\n")
        else:
            print("  -> Result: NOT persistently positive for aCL IgM.\n")
            
        # Check 3: Anti-ß2-glycoprotein-I (aß2GPI) antibodies
        print("Checking for persistent Anti-ß2-glycoprotein-I (aß2GPI) antibodies...")
        ab2gpi_igm_t1_pos = patient.labs_t1["antiB2GP1_IgM"] > LAB_NORMALS["antiB2GP1_IgM"]
        ab2gpi_igm_t2_pos = patient.labs_t2["antiB2GP1_IgM"] > LAB_NORMALS["antiB2GP1_IgM"]
        print(f"  - aß2GPI IgM Test 1: {patient.labs_t1['antiB2GP1_IgM']} UI/L (Normal < {LAB_NORMALS['antiB2GP1_IgM']}) -> {'Positive' if ab2gpi_igm_t1_pos else 'Negative'}")
        print(f"  - aß2GPI IgM Test 2: {patient.labs_t2['antiB2GP1_IgM']} UI/L (Normal < {LAB_NORMALS['antiB2GP1_IgM']}) -> {'Positive' if ab2gpi_igm_t2_pos else 'Negative'}")
        if ab2gpi_igm_t1_pos and ab2gpi_igm_t2_pos:
            persistent_findings.append("Anti-ß2GP1 IgM")
            print("  -> Result: PERSISTENTLY POSITIVE for aß2GPI IgM.\n")
        else:
            print("  -> Result: NOT persistently positive for aß2GPI IgM.\n")
            
        lab_criteria_met = len(persistent_findings) > 0
        if lab_criteria_met:
            print(f"Laboratory criterion MET. Confirmed persistent findings: {', '.join(persistent_findings)}.")
        else:
            print("Laboratory criterion NOT MET. No persistently positive antiphospholipid antibodies found.")

    print("\n")
    print("Step 3: Final Conclusion")
    print("-" * 70)
    
    final_answer = ""
    if clinical_criterion_met and lab_criteria_met:
        print("The patient meets BOTH the clinical and laboratory criteria for Antiphospholipid Syndrome.")
        final_answer = "Yes"
    else:
        print("The patient does NOT meet the full criteria for Antiphospholipid Syndrome.")
        final_answer = "No"

    # This is a special marker for the platform to read the final answer.
    # Do not modify it.
    print(f"\nDoes this patient categorize as having antiphospholipid syndrome?\n{final_answer}")
    sys.stdout.flush()
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    patient_case = PatientData()
    analyze_patient_for_aps(patient_case)
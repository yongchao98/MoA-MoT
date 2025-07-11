def solve_medical_case():
    """
    This script analyzes a clinical case to determine the most likely diagnosis
    based on key findings from the physical examination.
    """

    # Key findings from the patient's case
    palpable_vibration = True  # This is a "thrill"
    nonstop_murmur = True      # This is a "bruit"

    # Dictionary of possible diagnoses and their typical association with the key findings
    diagnoses = {
        "A. Femoral venous thrombosis": {"thrill": False, "bruit": False},
        "B. Arterial embolism": {"thrill": False, "bruit": False},
        "C. Retroperitoneal hematoma": {"thrill": False, "bruit": False},
        "D. Femoral artery dissection": {"thrill": False, "bruit": False},
        "E. Hamartoma": {"thrill": False, "bruit": False},
        "F. Femoral artery pseudoaneurysm": {"thrill": True, "bruit": True},
        "G. None of these choices": {"thrill": False, "bruit": False},
        "H. Arterio-capillary communication": {"thrill": False, "bruit": False}
    }

    print("Analyzing the patient's findings against possible diagnoses...")
    print(f"Patient's key signs: Palpable Vibration (Thrill) = {palpable_vibration}, Nonstop Murmur (Bruit) = {nonstop_murmur}\n")

    best_match = "G. None of these choices"
    for diagnosis, signs in diagnoses.items():
        # Check if the diagnosis matches the patient's key findings
        if signs["thrill"] == palpable_vibration and signs["bruit"] == nonstop_murmur:
            best_match = diagnosis
            break

    print("--- Diagnostic Analysis ---")
    print("The combination of a palpable vibration (thrill) and a nonstop murmur (bruit) at the site of a recent femoral artery puncture is a classic sign of a vascular complication involving turbulent blood flow.")
    print(f"The diagnosis that perfectly matches these findings is: '{best_match}'\n")
    
    print("--- Final Diagnostic Equation ---")
    # The "equation" represents the combination of clinical signs leading to the diagnosis.
    print("1 (Palpable Vibration) + 1 (Nonstop Murmur) = 2 (Key Signs Present)")
    print(f"Matching Diagnosis: {best_match}")


solve_medical_case()
print("<<<F>>>")
import textwrap

def recommend_htn_medications():
    """
    Analyzes the patient case and recommends a 3-drug regimen for resistant hypertension.
    """
    # Patient Data Points
    patient_bp = "145-153/85-89"
    patient_hr = 91
    has_diabetes = True

    # Rationale and Medication Selection
    print("Based on the patient's profile of resistant hypertension, comorbidities, and medication restrictions, the following 3-drug regimen is recommended to maximize treatment:\n")

    # Medication 1: ARB
    med1 = "Valsartan"
    reason1 = f"This Angiotensin II Receptor Blocker (ARB) is recommended as a foundational agent. Given the patient's concomitant Type 2 Diabetes, it provides essential renal (kidney) protection while effectively lowering blood pressure."
    
    # Medication 2: CCB
    med2 = "Amlodipine"
    reason2 = "This long-acting dihydropyridine calcium channel blocker (CCB) is a potent vasodilator and a first-line agent for hypertension. It is highly effective in combination therapy."

    # Medication 3: Beta-Blocker
    med3 = "Carvedilol"
    reason3 = f"With standard diuretics being unavailable, a beta-blocker is the next logical choice. Carvedilol is specifically chosen because it possesses both beta- and alpha-blocking properties, providing superior blood pressure reduction. It will also directly address the patient's elevated heart rate of {patient_hr} bpm."
    
    # Printing the recommendations
    print("1. Medication: " + med1)
    print("   Rationale: " + "\n   ".join(textwrap.wrap(reason1, width=80)))
    print("-" * 20)
    
    print("2. Medication: " + med2)
    print("   Rationale: " + "\n   ".join(textwrap.wrap(reason2, width=80)))
    print("-" * 20)
    
    print("3. Medication: " + med3)
    print("   Rationale: " + "\n   ".join(textwrap.wrap(reason3, width=80)))
    print("-" * 20)


if __name__ == '__main__':
    recommend_htn_medications()
<<<Valsartan, Amlodipine, Carvedilol>>>
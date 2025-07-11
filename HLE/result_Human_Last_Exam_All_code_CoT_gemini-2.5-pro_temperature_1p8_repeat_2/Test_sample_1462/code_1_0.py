def analyze_newborn_vitals(pounds, ounces):
    """
    Analyzes the newborn's weight from the clinical vignette.
    Converts weight to a single unit (kg) and checks for macrosomia.
    """
    # Conversion factors
    OUNCES_IN_A_POUND = 16.0
    KG_PER_POUND = 0.453592

    # --- Weight Conversion Equation ---
    # Weight_kg = (pounds + (ounces / OUNCES_IN_A_POUND)) * KG_PER_POUND
    
    # Calculate total weight in kg
    total_pounds = pounds + (ounces / OUNCES_IN_A_POUND)
    total_kg = total_pounds * KG_PER_POUND

    print("Analyzing patient's weight...")
    print(f"Given weight: {pounds} lb {ounces} oz")
    
    # Per the instructions, printing each number in the final equation
    print("\n--- Weight Conversion Calculation ---")
    print(f"Pounds component: {pounds}")
    print(f"Ounces component: {ounces}")
    print(f"Conversion factor (ounces per pound): {OUNCES_IN_A_POUND}")
    print(f"Conversion factor (kg per pound): {KG_PER_POUND}")
    print(f"Final calculated weight: {total_kg:.2f} kg")

    # Clinical context
    if total_kg > 4.0:
        print("\nClinical finding: Weight indicates macrosomia (> 4.0 kg).")
        print("This is often associated with maternal diabetes, a risk factor for various congenital anomalies.")
    else:
        print("\nClinical finding: Weight is within the normal range.")

# Vitals from the problem description
patient_pounds = 12
patient_ounces = 1

analyze_newborn_vitals(patient_pounds, patient_ounces)

def analyze_lupus_nephritis_labs():
    """
    Analyzes hypothetical lab results for a patient with a suspected
    lupus nephritis flare leading to rapid renal decline.
    """
    # Define normal lab value ranges for reference
    c3_normal_range_mg_dl = (80, 160)
    c4_normal_range_mg_dl = (15, 45)
    creatinine_normal_range_mg_dl = (0.6, 1.2)

    # Hypothetical lab results for the patient during acute renal decline
    patient_creatinine = 9.2  # mg/dL (Very high, indicating renal failure)
    patient_c3 = 40          # mg/dL (Very low)
    patient_c4 = 5           # mg/dL (Very low)

    print("--- Patient Lab Report ---")
    print(f"Patient's Serum Creatinine: {patient_creatinine} mg/dL")
    print(f"Normal Creatinine Range: {creatinine_normal_range_mg_dl[0]}-{creatinine_normal_range_mg_dl[1]} mg/dL")
    print("Interpretation: The high creatinine confirms severe kidney failure, but it is a result of the damage, not the cause.\n")

    print(f"Patient's Complement C3 Level: {patient_c3} mg/dL")
    print(f"Normal C3 Range: {c3_normal_range_mg_dl[0]}-{c3_normal_range_mg_dl[1]} mg/dL")
    is_c3_low = patient_c3 < c3_normal_range_mg_dl[0]
    print(f"Result: C3 level is critically low: {is_c3_low}\n")

    print(f"Patient's Complement C4 Level: {patient_c4} mg/dL")
    print(f"Normal C4 Range: {c4_normal_range_mg_dl[0]}-{c4_normal_range_mg_dl[1]} mg/dL")
    is_c4_low = patient_c4 < c4_normal_range_mg_dl[0]
    print(f"Result: C4 level is critically low: {is_c4_low}\n")

    print("--- Final Conclusion ---")
    if is_c3_low and is_c4_low:
        print("The critically low levels of complement C3 and C4 are the best indicators of the underlying cause.")
        print("They demonstrate massive consumption of complement proteins by immune complexes deposited in the kidneys.")
        print("This process, known as active lupus nephritis, is directly responsible for the rapid decline in renal function.")

# Execute the analysis
analyze_lupus_nephritis_labs()
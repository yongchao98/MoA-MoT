def diagnose_turner_phenocopy():
    """
    Analyzes clinical findings to determine the likely molecular abnormality
    in a patient with Turner-like features but a normal karyotype.
    """

    # --- Patient Data from Vignette ---
    patient_age = 18
    height_percentile = "10th - 15th"
    features = {
        "Amenorrhea and Infertility": True,
        "Ovarian Dysgenesis": True,
        "Short Stature": True,
        "Cardiovascular Symptoms (Fatigue, SOB, BP)": True,
        "Normal Karyotype (46,XX)": True
    }

    # --- Diagnostic Logic ---
    print("Step 1: Summarize key clinical and genetic findings.")
    print(f" - Patient Age: {patient_age} years")
    print(f" - Stature: Persistently in the {height_percentile} percentile.")
    for feature, present in features.items():
        if present:
            print(f" - Finding: {feature}")

    print("\nStep 2: Evaluate the combination of findings.")
    if features["Short Stature"] and features["Ovarian Dysgenesis"] and features["Cardiovascular Symptoms (Fatigue, SOB, BP)"]:
        print(" -> The constellation of short stature, ovarian failure, and cardiovascular symptoms is highly suggestive of Turner Syndrome.")
    
    print("\nStep 3: Incorporate the karyotype result.")
    if features["Normal Karyotype (46,XX)"]:
        print(" -> However, a normal karyotype rules out classic Turner Syndrome (which is 45,X0).")
        print(" -> This indicates a 'Turner Phenocopy' - a condition that mimics Turner Syndrome.")
        
    print("\nStep 4: Conclude the most likely molecular diagnosis.")
    print(" -> The most common cause of a Turner phenocopy is Noonan Syndrome.")
    print(" -> Noonan Syndrome is most frequently caused by a mutation in a specific gene.")
    
    final_gene = "PTPN11"
    print("\n-------------------------------------------------")
    print(f"Final Conclusion:")
    print(f"The likely molecular abnormality is a mutation in the {final_gene} gene.")
    print("-------------------------------------------------")

# Run the diagnostic analysis
diagnose_turner_phenocopy()
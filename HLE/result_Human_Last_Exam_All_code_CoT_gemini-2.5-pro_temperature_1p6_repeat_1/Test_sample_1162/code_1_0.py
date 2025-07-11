def diagnose_patient():
    """
    Analyzes clinical findings to determine the most likely diagnosis.
    """
    # Initialize the potential diagnoses with a likelihood score of 0.
    diagnoses = {
        "A. Germ cell tumor": 0,
        "B. Astrocytoma": 0,
        "C. Neuroblastoma": 0,
        "D. Nephroblastoma": 0,
        "E. Ewing sarcoma": 0
    }

    print("Analyzing the patient's clinical findings...\n")

    # Finding 1: Aniridia
    # This is a key finding, strongly associated with WAGR syndrome (Wilms tumor, Aniridia, Genitourinary anomalies, Range of developmental delays)
    print("Symptom: Aniridia")
    print("Logic: Strongly points to WAGR syndrome, which includes Nephroblastoma (Wilms tumor).")
    diagnoses["D. Nephroblastoma"] += 5
    print("Score added: Nephroblastoma + 5\n")


    # Finding 2: Pelvic Mass
    # A Nephroblastoma (Wilms tumor) is a kidney tumor that often presents as an abdominal or pelvic mass.
    print("Symptom: Pelvic Mass")
    print("Logic: Consistent with a large Nephroblastoma.")
    diagnoses["D. Nephroblastoma"] += 2
    diagnoses["A. Germ cell tumor"] += 1
    diagnoses["C. Neuroblastoma"] += 1
    print("Score added: Nephroblastoma + 2, Germ cell tumor + 1, Neuroblastoma + 1\n")


    # Finding 3: Hypertension (Elevated Blood Pressure)
    # Nephroblastoma can cause hypertension by secreting renin.
    print("Symptom: Hypertension")
    print("Logic: A known complication of Nephroblastoma (renin) and Neuroblastoma (catecholamines).")
    diagnoses["D. Nephroblastoma"] += 2
    diagnoses["C. Neuroblastoma"] += 1
    print("Score added: Nephroblastoma + 2, Neuroblastoma + 1\n")

    # Finding 4: Developmental Delay (Delayed Speech)
    # This fits the "R" (Range of developmental delays) in WAGR syndrome.
    print("Symptom: Developmental Delay")
    print("Logic: Associated with WAGR syndrome.")
    diagnoses["D. Nephroblastoma"] += 2
    print("Score added: Nephroblastoma + 2\n")

    # Calculate and display the final scores
    print("--- Final Likelihood Scores ---")
    most_likely_diagnosis = ""
    highest_score = -1

    # Using an equation-style format for the final tally
    print("Total for D. Nephroblastoma = 5 (Aniridia) + 2 (Mass) + 2 (Hypertension) + 2 (Delay)")
    
    total_nephroblastoma = 5 + 2 + 2 + 2
    diagnoses["D. Nephroblastoma"] = total_nephroblastoma
    
    for diagnosis, score in diagnoses.items():
        print(f"{diagnosis}: {score}")
        if score > highest_score:
            highest_score = score
            most_likely_diagnosis = diagnosis
    
    print("\n--- Conclusion ---")
    print(f"The most likely diagnosis based on the scoring is: {most_likely_diagnosis}")


if __name__ == "__main__":
    diagnose_patient()
<<<D>>>
import sys

def solve_diagnosis():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring potential options against key features of the case.
    """

    # Define the potential diagnoses
    diagnoses_options = {
        "A. Malignant Intertrigo": {"score": 0, "calc_list": []},
        "B. Allergic contact dermatitis": {"score": 0, "calc_list": []},
        "C. Hidradenitis Suppurativa": {"score": 0, "calc_list": []},
        "D. Atopic dermatitis": {"score": 0, "calc_list": []},
        "E. Psoriasis": {"score": 0, "calc_list": []}
    }

    print("Step 1: Analyzing clinical features and assigning likelihood points...\n")

    # --- Feature 1: Location (Axillary, Inframammary, Inguinal folds) ---
    # This distribution is classic for HS and Inverse Psoriasis.
    print("Feature: Location (Axillary, Inframammary, Inguinal folds)")
    diagnoses_options["C. Hidradenitis Suppurativa"]["calc_list"].append(3) # Classic triad
    diagnoses_options["E. Psoriasis"]["calc_list"].append(2) # Common for inverse psoriasis
    diagnoses_options["B. Allergic contact dermatitis"]["calc_list"].append(1) # Less specific
    diagnoses_options["D. Atopic dermatitis"]["calc_list"].append(1) # Less specific
    diagnoses_options["A. Malignant Intertrigo"]["calc_list"].append(1) # Possible but less classic distribution
    print("  - Hidradenitis Suppurativa: +3 (Classic)")
    print("  - Psoriasis: +2 (Common)")
    print("  - Others: +1 (Possible)\n")

    # --- Feature 2: Key Lesion (Purulent Nodules) ---
    # This is a hallmark of Hidradenitis Suppurativa.
    print("Feature: Lesion Type (Purulent Nodules)")
    diagnoses_options["C. Hidradenitis Suppurativa"]["calc_list"].append(3) # Hallmark
    diagnoses_options["E. Psoriasis"]["calc_list"].append(0) # Not a feature
    diagnoses_options["B. Allergic contact dermatitis"]["calc_list"].append(0) # Not a feature
    diagnoses_options["D. Atopic dermatitis"]["calc_list"].append(0) # Not a feature
    diagnoses_options["A. Malignant Intertrigo"]["calc_list"].append(0) # Not a feature
    print("  - Hidradenitis Suppurativa: +3 (Hallmark feature)")
    print("  - Others: +0 (Not a primary feature)\n")

    # --- Feature 3: Key Lesion (Large Bullae/Cysts) ---
    # Large abscesses/cysts in HS can be described as bullae.
    print("Feature: Lesion Type (Large Bullae/Cysts)")
    diagnoses_options["C. Hidradenitis Suppurativa"]["calc_list"].append(2) # Consistent with large abscesses
    diagnoses_options["B. Allergic contact dermatitis"]["calc_list"].append(1) # Possible in severe cases
    diagnoses_options["A. Malignant Intertrigo"]["calc_list"].append(1) # Can occur
    diagnoses_options["E. Psoriasis"]["calc_list"].append(0) # Not a feature
    diagnoses_options["D. Atopic dermatitis"]["calc_list"].append(0) # Not a feature
    print("  - Hidradenitis Suppurativa: +2 (Consistent)")
    print("  - Allergic Dermatitis / Malignant Intertrigo: +1 (Possible)")
    print("  - Others: +0 (Not a feature)\n")

    # --- Feature 4: Risk Factors (Obesity BMI 39, Smoking) ---
    # These are major, well-established risk factors for HS.
    print("Feature: Risk Factors (Obesity & Smoking)")
    diagnoses_options["C. Hidradenitis Suppurativa"]["calc_list"].append(3) # Major risk factors
    diagnoses_options["E. Psoriasis"]["calc_list"].append(1) # Associated, but less strongly
    diagnoses_options["B. Allergic contact dermatitis"]["calc_list"].append(0) # Not associated
    diagnoses_options["D. Atopic dermatitis"]["calc_list"].append(0) # Not associated
    diagnoses_options["A. Malignant Intertrigo"]["calc_list"].append(0) # Not associated
    print("  - Hidradenitis Suppurativa: +3 (Strongly associated)")
    print("  - Psoriasis: +1 (Weakly associated)")
    print("  - Others: +0 (Not specific risk factors)\n")

    # --- Calculate Final Scores ---
    print("---------------------------------------------------------")
    print("Step 2: Calculating final likelihood scores...")
    print("---------------------------------------------------------")

    most_likely_diagnosis = ""
    highest_score = -1

    for diagnosis, data in diagnoses_options.items():
        # Ensure list is not empty before summing
        if data["calc_list"]:
            total_score = sum(data["calc_list"])
            data["score"] = total_score

            # Format the equation string, handling the single case for clarity
            equation = " + ".join(map(str, data["calc_list"]))
            print(f"Score for {diagnosis}:")
            # This line fulfills the "output each number in the final equation" requirement
            print(f"  {equation} = {total_score}")

            if total_score > highest_score:
                highest_score = total_score
                most_likely_diagnosis = diagnosis
        else:
            print(f"Score for {diagnosis}: Not enough features to score.")


    print("\n---------------------------------------------------------")
    print("Step 3: Conclusion")
    print("---------------------------------------------------------")
    print(f"The most likely diagnosis based on this scoring system is '{most_likely_diagnosis}'.")
    print("\nReasoning: The combination of lesions in the classic intertriginous sites (axilla, inframammary, inguinal), the presence of purulent nodules (a hallmark), and major risk factors (obesity and smoking) makes Hidradenitis Suppurativa the strongest diagnosis.")

if __name__ == "__main__":
    solve_diagnosis()
<<<C>>>
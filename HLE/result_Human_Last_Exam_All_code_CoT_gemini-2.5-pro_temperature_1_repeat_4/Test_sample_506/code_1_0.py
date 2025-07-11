def recommend_htn_medications():
    """
    This function recommends and prints a 3-drug regimen for a patient with
    resistant hypertension, based on their clinical profile and contraindications.
    """
    # Recommended medications based on the analysis
    medication_1 = "Lisinopril"
    medication_2 = "Amlodipine"
    medication_3 = "Spironolactone"

    # Create a list of the medications
    recommended_regimen = [medication_1, medication_2, medication_3]

    # Print the recommendations
    print("Based on the patient's profile and medication restrictions, here are three recommended medications to maximize hypertension treatment:")
    for i, med in enumerate(recommended_regimen, 1):
        print(f"{i}. {med}")

# Execute the function to get the recommendations
recommend_htn_medications()
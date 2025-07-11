import sys

def solve_medical_case():
    """
    Analyzes a clinical case to provide a diagnosis, presenting the reasoning
    and incorporating patient data into an illustrative equation as requested.
    """
    # Step 1: Extract and list numerical data from the patient's case.
    age = 64
    para = 4
    gravida = 1
    bmi = 39
    cigarettes_daily = 3  # Using the upper end of the '2-3' range
    smoking_years = 15

    # Step 2: Print the analysis of the clinical findings.
    print("Analyzing the clinical case based on key features:")
    print("1. Primary Lesion: The presence of 'purulent nodules in the inguinal regions' is highly specific. This points to a deep-seated inflammatory process involving abscess formation.")
    print("2. Location: Lesions are in the axillary, inframammary, and inguinal folds. This pattern is characteristic of diseases affecting apocrine gland-bearing skin.")
    print("3. Patient Risk Factors: Two major risk factors are present:")
    print(f"   - Obesity (BMI of {bmi}): A well-established risk factor for Hidradenitis Suppurativa (HS).")
    print(f"   - Smoking ({cigarettes_daily} cigarettes/day for {smoking_years} years): The single most significant lifestyle risk factor associated with HS.")

    # Step 3: Evaluate the most likely diagnosis.
    print("\nEvaluating the differential diagnoses:")
    print("- Allergic/Atopic Dermatitis & Psoriasis: These conditions do not typically feature purulent nodules and deep abscesses as a primary sign.")
    print("- Hidradenitis Suppurativa (C): This chronic inflammatory disease perfectly matches the clinical picture: recurrent purulent nodules, abscesses, and sinus tracts in intertriginous areas, strongly associated with obesity and smoking.")

    # Step 4: Formulate the conclusion with the requested numerical 'equation'.
    print("\nConclusion:")
    print("The combination of symptoms and risk factors strongly supports the diagnosis of Hidradenitis Suppurativa.")
    print("\nAn illustrative equation showing the contribution of the patient's numerical data to the diagnosis is:")
    print(f"Diagnostic Confidence = f(Purulent Nodules) + f(BMI={bmi}) + f(Smoking={cigarettes_daily}cigs/day x {smoking_years}yrs) + f(Age={age}) + f(Para={para}, Gravida={gravida})")

# Execute the function to print the analysis.
solve_medical_case()
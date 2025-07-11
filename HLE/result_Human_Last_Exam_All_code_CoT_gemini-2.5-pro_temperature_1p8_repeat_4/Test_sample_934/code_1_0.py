def analyze_clinical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis
    by scoring risk factors and clinical findings.
    """
    # Patient demographic and social history values
    age = 53
    smoking_packs_per_day = 2
    smoking_years = 20
    
    # Calculate pack-years, a key metric for smoking history
    pack_years = smoking_packs_per_day * smoking_years

    # Assign scores to key findings based on their association with Esophageal Squamous Cell Carcinoma (SCC)
    # Higher scores indicate a stronger association.
    risk_scoring = {
        "Heavy Smoking (>40 pack-years)": 3,
        "Alcohol Use Disorder": 2,
        "Odynophagia (Painful Swallowing)": 2,
        "Imaging (Wall Thickening/Lumen Narrowing)": 2,
        "Age > 50": 1,
        "Endoscopy (Absence of ulcers/plaques)": 1 # This finding argues against other diagnoses like infection or GERD
    }
    
    # Prepare the components for the final output equation
    scores = list(risk_scoring.values())
    total_score = sum(scores)
    
    # Construct the equation string for clear output
    equation_string = " + ".join(map(str, scores))

    # Print the detailed analysis
    print("Analysis of Clinical Case:")
    print("==========================")
    print(f"Patient is a {age}-year-old woman with significant risk factors.")
    print(f"Smoking History: {pack_years} pack-years, a major risk factor for SCC.")
    print("\nA scoring model helps quantify the risk for Esophageal Squamous Cell Carcinoma:")
    for factor, score in risk_scoring.items():
        print(f"- {factor}: Score = {score}")

    print("\nCalculating the total risk score using the individual values:")
    print(f"Risk Score Equation: {equation_string} = {total_score}")

    print("\nConclusion:")
    print("The patient's combination of heavy, long-term smoking and alcohol use are the most significant risk factors for Esophageal Squamous Cell Carcinoma.")
    print("Her symptoms (pain with swallowing) and imaging results (wall thickening) strongly suggest an esophageal mass.")
    print("Critically, the endoscopy was negative for surface-level signs like ulcers or plaques, making infectious causes (Herpes, Strep) or GERD less likely. This points towards an infiltrating tumor growing within the esophageal wall.")
    print("\nGiven the overwhelming evidence from her risk profile and clinical presentation, the most likely diagnosis is Esophageal Squamous Cell Carcinoma.")

analyze_clinical_case()
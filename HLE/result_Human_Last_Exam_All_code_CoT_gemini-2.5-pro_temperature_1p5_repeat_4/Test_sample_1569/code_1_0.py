def solve_medical_case():
    """
    Analyzes a clinical vignette to determine the most likely diagnosis.
    This script uses a weighted scoring system to represent the diagnostic process.
    """

    # Assign diagnostic weights to the key findings in the case.
    # The most specific finding receives the highest weight.
    symptoms_weight = 1
    geographic_risk_weight = 1
    disseminated_disease_weight = 2
    positive_lyme_serology_weight = 10 # This is the most definitive finding.

    # The positive "Lyme serology" directly points to the causative agent of Lyme disease.
    causative_agent = "Borrelia burgdorferi"

    # The final diagnostic "equation" is the sum of these weights.
    total_score = symptoms_weight + geographic_risk_weight + disseminated_disease_weight + positive_lyme_serology_weight

    print("Analyzing the clinical case based on a weighted diagnostic model:")
    print(f"- Nonspecific symptoms (fever, headaches): {symptoms_weight} point")
    print(f"- Geographic risk (Oklahoma camping): {geographic_risk_weight} point")
    print(f"- Disseminated signs (disorientation, murmur): {disseminated_disease_weight} points")
    print(f"- Positive Lyme serology (elevated IgM): {positive_lyme_serology_weight} points")
    print("-" * 30)

    # Print the equation as requested
    print("The diagnostic score equation is:")
    print(f"{symptoms_weight} + {geographic_risk_weight} + {disseminated_disease_weight} + {positive_lyme_serology_weight} = {total_score}")
    print("-" * 30)

    print("Conclusion:")
    print("The lab result explicitly states a positive 'Lyme serology titer' for an acute infection (elevated IgM).")
    print(f"The organism that causes Lyme disease is {causative_agent}.")
    print(f"Therefore, the positive titer is for {causative_agent}.")
    print("\nAnswer Choice C corresponds to Borrelia burgdorferi.")

solve_medical_case()
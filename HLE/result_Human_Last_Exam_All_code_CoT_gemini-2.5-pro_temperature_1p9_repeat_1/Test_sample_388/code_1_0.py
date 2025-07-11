def generate_counseling_recommendation():
    """
    Analyzes a patient's medications to identify interactions and generate a
    pharmacist's counseling recommendation.
    """
    # Medications List
    otc_med = "Excedrin (contains Aspirin, an NSAID)"
    prescription_1 = "Atorvastatin 20mg"
    prescription_2 = "Junel Fe 1.5/30mg"
    prescription_3 = "Fluoxetine 20mg (an SSRI)"

    # Identifying the specific drugs involved in the interaction
    interacting_drug_1 = "Fluoxetine 20mg"
    interacting_drug_2 = "Aspirin (in Excedrin)"

    # Formulating the recommendation
    print("Pharmacist Counseling Recommendation")
    print("-" * 40)
    print("Patient: Allison")
    print(f"Medications on file: {prescription_1}, {prescription_2}, {prescription_3}")
    print(f"Patient also reports taking: {otc_med}\n")

    print("--- POTENTIAL INTERACTION IDENTIFIED ---")
    print(f"Interaction between: {interacting_drug_1} and {interacting_drug_2}.")
    print("\nRisk Analysis:")
    print("Both SSRIs (like Fluoxetine) and NSAIDs (like Aspirin) can independently increase the risk of gastrointestinal bleeding.")
    print("When taken together, this risk is significantly increased.\n")

    print("--- RECOMMENDED ACTION ---")
    print("The pharmacist should provide the following counseling advice:")
    print("\n\t'Hello Allison, I see you are taking Fluoxetine 20mg and also using Excedrin.")
    print("\tI want to let you know that taking Fluoxetine with the aspirin in Excedrin can increase")
    print("\tyour risk of stomach bleeding. For future headaches or pain, a safer alternative")
    print("\twould be a product with just acetaminophen, such as Tylenol.'")
    print("\n\t'If you ever notice symptoms like black, tarry stools, stomach pain that doesn't go away,")
    print("\tor vomiting blood, please contact your doctor right away.'")
    print("-" * 40)

# Execute the function to print the recommendation
generate_counseling_recommendation()
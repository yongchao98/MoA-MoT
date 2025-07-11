def pharmacist_counseling():
    """
    Provides a counseling recommendation based on Allison's medications.
    """
    medication_1 = "Fluoxetine 20mg"
    medication_2 = "Excedrin"
    component_of_med_2 = "Aspirin"
    
    drug_class_1 = "SSRI (Selective Serotonin Reuptake Inhibitor)"
    drug_class_2 = "NSAID (Nonsteroidal Anti-inflammatory Drug)"
    
    interaction = "increased risk of gastrointestinal bleeding"
    recommendation = "use an alternative pain reliever like acetaminophen (Tylenol) for headaches"
    
    print("Pharmacist Counseling Recommendation:")
    print("-" * 35)
    print(f"There is a potential interaction between your new prescription, {medication_1}, and the {medication_2} you took for your headache.")
    print(f"Explanation: {medication_1} is an {drug_class_1}, and {medication_2} contains {component_of_med_2}, which is an {drug_class_2}.")
    print(f"When taken together, these two types of drugs can significantly increase the {interaction}.")
    print(f"\nRecommendation: We recommend you {recommendation} in the future while you are taking {medication_1}.")
    print("Please discuss this with your doctor to determine the best and safest option for you.")

pharmacist_counseling()
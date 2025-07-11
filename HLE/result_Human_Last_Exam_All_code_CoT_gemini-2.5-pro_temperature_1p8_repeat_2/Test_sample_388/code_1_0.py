def generate_counseling_recommendation():
    """
    Analyzes Allison's medications and generates a pharmacist's counseling recommendation.
    """
    # Define medications and their relevant classes
    fluoxetine_mg = 20
    excedrin_components = "Aspirin, Acetaminophen, and Caffeine"
    fluoxetine_class = "SSRI (Selective Serotonin Reuptake Inhibitor)"
    aspirin_class = "NSAID (Nonsteroidal Anti-inflammatory Drug)"

    # Print the recommendation step-by-step
    print("Pharmacist Counseling Recommendation")
    print("=" * 40)
    print("A review of your medications shows a potential interaction that is important to discuss.\n")

    print(f"1. Prescription Medication: Fluoxetine {fluoxetine_mg}mg")
    print(f"   - Class: {fluoxetine_class}\n")
    print(f"2. Over-the-Counter Medication: Excedrin")
    print(f"   - Active Ingredient of Concern: Aspirin, which is an {aspirin_class}\n")

    print("Interaction Concern:")
    print(f"When Fluoxetine (an {fluoxetine_class}) and Aspirin (an {aspirin_class}) are taken together,")
    print("they can significantly increase the risk of bleeding, especially in the stomach and digestive tract.\n")

    print("Recommendation:")
    print("a. For future headaches or pain, it is recommended to choose a product that only contains")
    print("   acetaminophen (such as Tylenol) and avoid Excedrin or other NSAIDs (like ibuprofen or naproxen).\n")
    print("b. Please be sure to monitor for any signs of unusual bleeding, such as easy bruising,")
    print("   prolonged bleeding from cuts, nosebleeds, or dark, tarry stools.")
    print("   Contact your doctor if you notice any of these symptoms.\n")

# Run the function to display the counseling recommendation
generate_counseling_recommendation()
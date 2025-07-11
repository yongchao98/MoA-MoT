def provide_counseling():
    """
    This function formulates and prints a pharmacist's counseling recommendation
    based on a potential drug interaction.
    """
    # Medications involved in the interaction
    prescription_drug = "Fluoxetine"
    otc_drug = "Excedrin"
    interacting_component = "aspirin"
    drug_class = "NSAID (Nonsteroidal Anti-Inflammatory Drug)"

    # Details of the interaction and recommendation
    risk = "an increased risk of stomach bleeding"
    safer_alternative = "acetaminophen (Tylenol)"

    # Print the counseling point in a clear, conversational format
    print("--- Pharmacist Counseling Recommendation ---")
    print(f"\nHello Allison, while reviewing your medications, I see you are taking {prescription_drug}.")
    print(f"You mentioned you took {otc_drug} for your headache. It's important for you to be aware of a potential interaction.")
    print(f"\nThe {interacting_component} in {otc_drug}, which is an {drug_class}, can interact with {prescription_drug}.")
    print(f"\nTaking these two together can result in {risk}.")
    print("\nMy recommendation is to use a safer option for pain relief in the future.")
    print(f"A better choice would be a product with just {safer_alternative}, which does not have this interaction.")
    print("\nPlease let your doctor know about your headaches and always feel free to ask me before starting any new over-the-counter medications.")
    print("\n--- End of Recommendation ---")

provide_counseling()
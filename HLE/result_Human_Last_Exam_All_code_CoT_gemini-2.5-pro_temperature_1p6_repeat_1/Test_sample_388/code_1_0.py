import textwrap

def generate_counseling_recommendation():
    """
    This function generates and prints a pharmacy counseling recommendation
    based on Allison's medications and symptoms.
    """
    
    title = "Pharmacist Counseling Recommendation"
    print(title)
    print("=" * len(title))

    print("\nPatient Medications Identified:")
    print("- Prescription: Atorvastatin 20mg")
    print("- Prescription: Junel Fe 1.5/30mg")
    print("- Prescription: Fluoxetine 20mg")
    print("- Over-the-Counter (OTC): Excedrin (contains Aspirin, an NSAID)")

    print("\n--- CRITICAL COUNSELING POINT: DRUG INTERACTION ---")
    
    interaction_header = "\nInteraction Detected:"
    interaction_detail = "Fluoxetine (an SSRI antidepressant) and Aspirin (an NSAID found in Excedrin)."
    print(interaction_header)
    print(interaction_detail)

    risk_header = "\nAssociated Risk:"
    risk_detail = "Both Fluoxetine and Aspirin increase the risk of bleeding, especially in the stomach and gastrointestinal tract. Using them together significantly magnifies this risk."
    print(risk_header)
    for line in textwrap.wrap(risk_detail, width=70):
        print(line)

    recommendation_header = "\nPrimary Recommendation:"
    recommendation_1 = "1. For future headaches, it is strongly recommended to use a safer alternative that does not contain an NSAID. Acetaminophen (the main ingredient in Tylenol) taken alone is a much safer choice."
    recommendation_2 = "2. Discontinue the use of Excedrin while on Fluoxetine therapy to avoid the risk of bleeding."
    recommendation_3 = "3. Be aware of the signs of gastrointestinal bleeding, such as black or tarry stools, severe stomach pain, or vomiting blood or material that looks like coffee grounds, and seek medical attention if they occur."

    print(recommendation_header)
    for line in textwrap.wrap(recommendation_1, width=70):
        print(line)
    for line in textwrap.wrap(recommendation_2, width=70):
        print(line)
    for line in textwrap.wrap(recommendation_3, width=70):
        print(line)


if __name__ == "__main__":
    generate_counseling_recommendation()
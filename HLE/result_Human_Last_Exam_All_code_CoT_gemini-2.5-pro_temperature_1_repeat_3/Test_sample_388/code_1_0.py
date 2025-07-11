import textwrap

def provide_counseling_recommendation():
    """
    Analyzes Allison's medications to identify the most critical interaction
    and prints a clear counseling recommendation for the pharmacist to make.
    """
    # Define the relevant medications from the scenario
    ssri_medication = "Fluoxetine 20mg"
    otc_pain_reliever = "Excedrin"
    interacting_ingredient = "Aspirin"
    
    # Formulate the counseling points based on the known interaction
    # between SSRIs (Fluoxetine) and NSAIDs (Aspirin).
    
    title = "Pharmacist Counseling Recommendation for Allison"
    
    introduction = (
        f"Based on your prescriptions and your recent use of {otc_pain_reliever}, "
        "there is one very important counseling point to discuss."
    )
    
    interaction_explanation = (
        f"Your prescription, {ssri_medication}, belongs to a class of drugs called SSRIs. The {otc_pain_reliever} "
        f"you took for your headache contains {interacting_ingredient}, which is an NSAID. "
        "When an SSRI and an NSAID are taken together, there is a significantly "
        "increased risk of bleeding, especially in the stomach."
    )
    
    recommendation = (
        "For future headaches or pain, it is much safer for you to choose a medication "
        "that contains only acetaminophen (such as Tylenol). Acetaminophen "
        "does not have this interaction and is a safer choice for you while you are taking Fluoxetine."
    )
    
    monitoring = (
        "You should watch for any signs of gastrointestinal bleeding, such as "
        "unusual stomach pain, black or tarry-looking stools, or vomiting anything that "
        "looks like coffee grounds. If you notice any of these, please contact your "
        "doctor right away."
    )
    
    # Print the final, formatted recommendation
    print(f"\n{title}")
    print("-" * len(title))
    
    print("\n[Interaction Alert]")
    print(textwrap.fill(interaction_explanation, width=80))
    
    print("\n[Recommendation]")
    print(textwrap.fill(recommendation, width=80))
    
    print("\n[Important Monitoring]")
    print(textwrap.fill(monitoring, width=80))

# Execute the function to print the recommendation
provide_counseling_recommendation()
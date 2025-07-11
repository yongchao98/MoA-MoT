def generate_counseling_recommendation():
    """
    Analyzes Allison's medications and generates a pharmacist's counseling recommendation.
    """

    # Define the medications Allison is taking
    fluoxetine_dose = 20
    junel_fe_dose = "1.5/30"
    atorvastatin_dose = 20
    otc_med = "Excedrin"
    otc_ingredient = "Aspirin"

    # Define the core counseling points
    interaction_med1 = f"Fluoxetine {fluoxetine_dose}mg"
    interaction_med2 = f"{otc_med} (which contains {otc_ingredient})"
    interaction_risk = "an increased risk of stomach bleeding"
    safer_alternative = "acetaminophen (the active ingredient in Tylenol)"

    symptom_med = f"Junel Fe {junel_fe_dose}mg"
    symptom_alert = "new or worsening headaches can be a serious side effect of oral contraceptives"

    # Construct and print the final recommendation text
    recommendation = (
        f"Based on your prescriptions, the main counseling recommendation involves the {otc_med} you took for your headache.\n"
        f"1. There is an interaction between your {interaction_med1} and the {interaction_med2}. When taken together, they can create {interaction_risk}. For future headaches, it would be safer to use a product containing just {safer_alternative}.\n"
        f"2. Also, because you are taking {symptom_med}, it's important to know that {symptom_alert}. If your headaches continue, get worse, or feel different from normal, please make sure to contact your doctor.\n"
        f"Your Atorvastatin {atorvastatin_dose}mg does not have a significant interaction here."
    )

    print(recommendation)

generate_counseling_recommendation()

# The final answer is the text generated above.
final_answer = (
    "Based on your prescriptions, the main counseling recommendation involves the Excedrin you took for your headache.\n"
    "1. There is an interaction between your Fluoxetine 20mg and the Excedrin (which contains Aspirin). When taken together, they can create an increased risk of stomach bleeding. For future headaches, it would be safer to use a product containing just acetaminophen (the active ingredient in Tylenol).\n"
    "2. Also, because you are taking Junel Fe 1.5/30mg, it's important to know that new or worsening headaches can be a serious side effect of oral contraceptives. If your headaches continue, get worse, or feel different from normal, please make sure to contact your doctor.\n"
    "Your Atorvastatin 20mg does not have a significant interaction here."
)
print(f"\n<<<{final_answer}>>>")
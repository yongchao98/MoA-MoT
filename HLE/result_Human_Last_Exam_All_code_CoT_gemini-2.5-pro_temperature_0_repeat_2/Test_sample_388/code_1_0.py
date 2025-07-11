def pharmacist_counseling():
    """
    Generates a counseling recommendation based on Allison's medications.
    """
    patient_name = "Allison"
    ssri_drug = "Fluoxetine"
    otc_headache_med = "Excedrin"
    interacting_ingredient = "Aspirin"
    safer_alternative = "acetaminophen (like Tylenol)"

    recommendation = (
        f"Hello {patient_name}, I see you're picking up your {ssri_drug}. "
        f"I also understand you recently took {otc_headache_med} for a headache.\n\n"
        f"It's important for you to know that taking {ssri_drug} together with {interacting_ingredient}, "
        f"which is an active ingredient in {otc_headache_med}, can increase your risk of bleeding, especially in the stomach.\n\n"
        f"For future headaches, it would be a safer choice to use a product that only contains {safer_alternative}.\n\n"
        f"Please be sure to watch for any signs of unusual bleeding or bruising, like nosebleeds or dark stools, and contact your doctor if you notice any."
    )

    print(recommendation)
    
    # The final answer is the full text of the recommendation.
    final_answer = f"<<<{recommendation}>>>"
    # The line below is for the final answer extraction and is not part of the counseling text.
    # print(f"\n\n{final_answer}")

if __name__ == "__main__":
    pharmacist_counseling()
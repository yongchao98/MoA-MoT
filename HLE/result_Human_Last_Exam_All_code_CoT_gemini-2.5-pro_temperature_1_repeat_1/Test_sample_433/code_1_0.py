def get_bud_for_single_dose_ampule():
    """
    This function provides the Beyond-Use Date (BUD) for a single-dose ampule
    based on United States Pharmacopeia (USP) Chapter <797> guidelines.
    """
    # According to USP <797>, a sterile environment is ISO Class 5 or cleaner.
    # The rule for a single-dose container (like an ampule) opened in such an environment is specific.
    bud_hours = 12
    container = "single-dose ampule"
    environment = "sterile environment (ISO Class 5 or better)"

    print(f"Guideline Source: USP General Chapter <797> Pharmaceutical Compounding—Sterile Preparations")
    print(f"Container Type: {container}")
    print(f"Compounding Environment: {environment}")
    print("-" * 50)
    print(f"From the time of puncture, the Beyond-Use Date (BUD) is {bud_hours} hours.")
    print("-" * 50)
    print("This means the medication must be administered or discarded within this timeframe.")

# Execute the function to display the information.
get_bud_for_single_dose_ampule()
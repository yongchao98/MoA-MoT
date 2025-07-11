def get_bud_for_single_dose_ampule():
    """
    This function explains the Beyond-Use Date (BUD) for a single-dose ampule
    according to USP <797> guidelines.
    """
    # According to USP <797>, the BUD for a single-dose container opened
    # in an ISO Class 5 environment is 12 hours.
    bud_in_hours = 12

    print("According to United States Pharmacopeia (USP) <797> guidelines for a single-dose container (e.g., an ampule) opened in a sterile, ISO Class 5 environment:")
    print("\nThe Beyond-Use Date (BUD) is determined from the time of puncture or opening.")

    # Display the final "equation" or rule
    print("\nThe final rule is:")
    print(f"Time from puncture + {bud_in_hours} hours = Expiration Time")
    print(f"\nTherefore, the BUD is {bud_in_hours} hours.")

get_bud_for_single_dose_ampule()
def get_bud_for_single_dose_ampule():
    """
    This function explains and provides the Beyond-Use Date (BUD) for a single-dose
    container like an ampule when opened in a sterile environment, based on
    USP <797> guidelines.
    """
    # According to USP <797>, a single-dose container is one intended for a single use
    # and typically lacks an antimicrobial preservative.
    # An ampule is a prime example.
    
    # The environment specified is "sterile," which corresponds to an
    # ISO Class 5 or cleaner environment in pharmacy practice.
    
    # The BUD for a single-dose container (vial, ampule, etc.) that has been
    # opened or punctured in an ISO Class 5 environment is 12 hours.
    bud_in_hours = 12
    
    print("Based on United States Pharmacopeia (USP) <797> guidelines for sterile compounding:")
    print("A single-dose container (such as an ampule or vial) that is punctured or opened in a sterile, ISO Class 5 environment must be used within the following timeframe.")
    print("\nFinal Answer Equation:")
    print(f"The Beyond-Use Date (BUD) is {bud_in_hours} hours after the time of initial puncture or opening.")
    print("\nNote: An ampule cannot be resealed after opening. It should be used immediately, and any unused portion must be discarded.")

get_bud_for_single_dose_ampule()
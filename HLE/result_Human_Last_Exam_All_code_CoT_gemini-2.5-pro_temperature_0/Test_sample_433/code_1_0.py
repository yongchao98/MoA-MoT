def get_bud_for_ampule():
    """
    Provides the Beyond-Use Date (BUD) for a single-dose ampule
    based on USP <797> guidelines.
    """
    # According to USP <797>, a sterile environment for compounding
    # is typically an ISO Class 5 environment.
    environment = "ISO Class 5 or better sterile environment"
    container_type = "single-dose container (e.g., an ampule)"
    
    # The BUD for a single-dose container opened in an ISO Class 5 environment.
    bud_hours = 12
    
    # Note: An ampule is a glass container that cannot be resealed after opening.
    # Therefore, the ampule itself must be discarded immediately after drawing up the medication.
    # The BUD applies to the medication that has been drawn from the ampule and maintained
    # in a sterile state (e.g., in a syringe).
    
    print(f"According to USP <797> guidelines:")
    print(f"The Beyond-Use Date (BUD) for the contents of a {container_type}")
    print(f"when opened in an {environment} is up to:")
    print(f"\n{bud_hours} hours from the time of puncture/opening.")
    print("\nImportant Note: An opened ampule cannot be stored and must be discarded immediately after use.")
    print("The BUD applies to the medication drawn from the ampule, not the opened container itself.")

get_bud_for_ampule()
def get_ampule_bud():
    """
    Calculates and explains the Beyond-Use Date (BUD) for a single-dose ampule.
    """
    container_type = "single-dose ampule"
    environment = "sterile environment"
    
    # According to pharmacy sterile compounding standards (e.g., USP <797>),
    # an ampule is a glass container that is broken open and cannot be resealed.
    # This design means that once opened, its contents are exposed to the environment.
    
    # To prevent contamination and ensure patient safety, the contents must be used immediately.
    # Any unused portion must be discarded and not stored for later use.
    # The term "immediately" can be numerically represented as a storage time of 0 hours.
    
    bud_in_hours = 0
    
    print(f"Container Type: {container_type}")
    print(f"Opening Condition: Punctured/opened in a {environment}")
    print("-" * 30)
    print("Conclusion on Beyond-Use Date (BUD):")
    print("The contents of an opened ampule must be used immediately.")
    print(f"The maximum time permitted for storage after opening is: {bud_in_hours} hours.")
    print("Any unused portion must be discarded immediately and cannot be saved.")

get_ampule_bud()
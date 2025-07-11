def get_ampule_bud():
    """
    Explains the Beyond-Use Date (BUD) for a single-dose ampule opened in a sterile environment.
    """
    container_type = "Single-Dose Ampule"
    key_principle = "Cannot be resealed after opening, compromising sterility."
    
    # According to USP <797> and ISMP guidelines, there is no storage time for an open ampule.
    bud_storage_time_hours = 0
    
    print(f"Container Type: {container_type}")
    print(f"Key Principle: {key_principle}")
    print("\n--- Official Guideline ---")
    print("Once a single-dose ampule is opened, its contents must be used immediately.")
    print("Any unused portion must be discarded and cannot be stored for later use.")
    print("\n--- Final Answer ---")
    print(f"The Beyond-Use Date (BUD) for storage after opening is: {bud_storage_time_hours} hours.")

get_ampule_bud()
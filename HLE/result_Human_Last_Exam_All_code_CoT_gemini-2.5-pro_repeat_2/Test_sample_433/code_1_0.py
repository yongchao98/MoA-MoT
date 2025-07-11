def get_ampule_bud():
    """
    This function determines and explains the Beyond-Use Date (BUD)
    for a single-dose ampule opened in a sterile environment.
    """
    container_type = "Single-dose ampule"
    environment = "Sterile environment (ISO Class 5 or cleaner)"
    
    # According to USP <797> guidelines, an ampule is an open system once broken.
    # It cannot be resealed, and its contents are for immediate use only.
    # There is no allowable storage time after opening.
    storage_time_hours = 0
    
    print(f"Container Type: {container_type}")
    print(f"Environment: {environment}")
    print("-" * 40)
    print("Guideline: Once an ampule is opened, it cannot be resealed.")
    print("Its contents must be used immediately, and any unused portion must be discarded.")
    print(f"The allowable storage time is {storage_time_hours} hours.")
    print("Therefore, the Beyond-Use Date (BUD) is 'Immediate Use'.")

get_ampule_bud()
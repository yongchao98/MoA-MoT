def get_ampule_bud():
    """
    This function explains the Beyond-Use Date (BUD) for a single-dose ampule
    according to USP <797> guidelines.
    """
    # According to USP <797>, an opened ampule cannot be stored.
    # It must be used immediately.
    container_type = "single-dose ampule"
    allowable_storage_time_hours = 0
    
    # The final rule is based on the inability to reseal an ampule.
    # Therefore, the equation for the BUD is: Time of Use = Time of Puncture + 0 hours.
    
    print(f"The Beyond-Use Date (BUD) for a {container_type} opened in a sterile environment is immediate.")
    print(f"It must be used immediately after opening, and any unused portion must be discarded.")
    print(f"The allowable storage time is {allowable_storage_time_hours} hours.")

get_ampule_bud()
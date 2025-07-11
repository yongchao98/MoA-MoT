def get_ampule_bud():
    """
    This function explains the Beyond-Use Date (BUD) for a single-dose ampule.
    """
    # According to USP <797> sterile compounding standards, an ampule is a single-dose container
    # that cannot be resealed after it is opened.
    
    # The rule for an opened ampule in a sterile environment is that it must be used immediately.
    # It cannot be stored for any period of time.
    bud_for_storage_in_hours = 0
    
    print("For a single-dose container ampule punctured in a sterile environment:")
    print("The contents must be used immediately.")
    print(f"The Beyond-Use Date (BUD) for storage is {bud_for_storage_in_hours} hours, as it cannot be stored.")
    print("After the initial dose is withdrawn, the ampule and any remaining contents must be discarded immediately.")

get_ampule_bud()
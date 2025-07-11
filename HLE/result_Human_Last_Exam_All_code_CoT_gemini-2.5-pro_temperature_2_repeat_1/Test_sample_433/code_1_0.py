def get_ampule_bud():
    """
    Determines and explains the Beyond-Use Date (BUD) for a single-dose ampule
    opened in a sterile environment, based on USP <797> guidelines.
    """

    # According to USP <797>, the allowable storage time for an opened ampule is zero.
    allowable_storage_time_hours = 0

    # --- Explanation ---
    print("According to United States Pharmacopeia (USP) General Chapter <797> guidelines:")
    print("1. An ampule is a single-dose container that cannot be resealed after being opened.")
    print("2. Because it is an open system, an opened ampule must not be stored for any period of time.")
    print("3. The contents must be withdrawn and used immediately.")

    # --- Equation Representation ---
    print("\nThis principle can be shown with the following equation:")
    print("Beyond-Use Date (BUD) = Time of Opening + Allowable Storage Time")

    print("\nFor a single-dose ampule:")
    # Here is the number from the final equation as requested
    print(f"Allowable Storage Time = {allowable_storage_time_hours} hours")
    
    # --- Final Conclusion ---
    print("\nTherefore, the Beyond-Use Date for the contents of the ampule is 'Immediate Use'.")

# Execute the function
get_ampule_bud()
def get_bud_information():
    """
    This function provides the Beyond-Use Date (BUD) information for a single-dose ampule
    and compares it with a single-dose vial, based on USP <797> guidelines for sterile environments.
    """

    # Information for a single-dose ampule
    container_ampule = "single-dose ampule"
    bud_ampule_hours = 0
    bud_ampule_rule = "Immediate use only"

    # Information for a single-dose vial (for comparison)
    container_vial = "single-dose vial"
    bud_vial_hours = 12
    bud_vial_rule = f"Up to {bud_vial_hours} hours after initial puncture"

    # --- Main Answer ---
    print(f"Analysis for: A {container_ampule} punctured in a sterile environment.")
    print("-" * 60)
    print("An ampule is a sealed glass container that is broken open and cannot be resealed.")
    print("Once opened, it is an open system exposed to the environment.")
    print("\nAccording to sterile compounding standards (USP <797>):")
    print(f"*   The BUD is: '{bud_ampule_rule}'.")
    print(f"*   Storage time allowed: {bud_ampule_hours} hours.")
    print("*   Any unused portion must be discarded immediately and cannot be stored for future use.")
    print("-" * 60)

    # --- Contextual Comparison ---
    print("\nFor comparison, the rule for a single-dose VIAL is different:")
    print(f"A {container_vial}, which has a rubber stopper, may be used for up to {bud_vial_hours} hours after the first puncture when maintained in a sterile environment.")

if __name__ == "__main__":
    get_bud_information()
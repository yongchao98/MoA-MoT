def get_bud_information():
    """
    Explains the Beyond-Use Date (BUD) for a single-dose ampule.
    """
    container_type = "single-dose ampule"
    action = "punctured/opened in a sterile environment"
    bud_ampule = "Immediate Use"

    # For comparison
    comparison_container = "single-dose vial"
    bud_vial_hours = 12

    print(f"Container Type: {container_type}")
    print(f"Action: {action}")
    print("-" * 30)
    print(f"The official Beyond-Use Date (BUD) is: {bud_ampule}")
    print("\nExplanation:")
    print("An ampule is a glass container that cannot be resealed once broken open.")
    print("Because it becomes an open system, the contents must be withdrawn and used immediately to prevent microbial contamination.")
    print("Any unused portion must be discarded and cannot be stored for any length of time.")
    
    print("\nFor Comparison:")
    print(f"A {comparison_container}, when punctured in the same sterile environment, has a BUD of up to {bud_vial_hours} hours.")
    print(f"This is because a vial's rubber stopper allows it to remain a closed system, unlike an open ampule.")

get_bud_information()
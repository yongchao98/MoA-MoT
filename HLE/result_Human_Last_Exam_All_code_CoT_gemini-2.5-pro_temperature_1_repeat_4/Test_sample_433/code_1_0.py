def get_bud_for_single_dose_container():
    """
    This function provides the Beyond-Use Date (BUD) for a single-dose
    container based on the environment in which it is punctured,
    according to USP <797> guidelines.
    """
    # BUD in hours for a single-dose container punctured in a sterile environment
    # (ISO Class 5 or cleaner).
    bud_in_sterile_env_hours = 12

    # For context, the BUD in hours for a container opened outside a sterile environment.
    bud_outside_sterile_env_hours = 1

    # The user's question specifies a sterile environment.
    print(f"The Beyond-Use Date (BUD) for a single-dose container ampule from the time of puncture in a sterile environment (e.g., an ISO Class 5 PEC) is up to {bud_in_sterile_env_hours} hours.")
    print("\nIt is important to note that if the same ampule were opened in an environment worse than ISO Class 5, the contents must be used within 1 hour or discarded.")

get_bud_for_single_dose_container()
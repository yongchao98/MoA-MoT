def get_bud_for_ampule():
    """
    This function provides the Beyond-Use Date (BUD) for a single-dose ampule
    punctured in a sterile environment, based on USP <797> guidelines.
    """
    # According to USP <797>, the BUD is 12 hours for a single-dose container
    # opened in an ISO Class 5 environment.
    bud_hours = 12
    number_of_doses = 1

    # The final "equation" is a descriptive sentence using the variables.
    print(f"The Beyond-Use Date (BUD) for a single ({number_of_doses}) dose container ampule, from the time of puncture in a sterile environment, is up to {bud_hours} hours.")
    print("\nThis assumes the ampule is opened and maintained within an ISO Class 5 or better environment.")

get_bud_for_ampule()
def get_bud_for_single_dose_ampule():
    """
    Calculates and prints the Beyond-Use Date (BUD) for a single-dose ampule
    opened in a sterile environment, based on USP <797> guidelines.
    """
    # Parameters based on the query
    container_type = "single-dose container ampule"
    environment = "sterile environment (e.g., ISO Class 5)"
    
    # The BUD according to USP <797> guidelines
    bud_in_hours = 12
    
    # The "equation" here is simply stating the resulting number for the BUD.
    # The final output will explicitly show the number 12.
    print(f"For a {container_type} opened in a {environment}, the Beyond-Use Date (BUD) is up to {bud_in_hours} hours after the time of puncture.")
    print(f"The final number for the BUD is: {bud_in_hours}")

if __name__ == "__main__":
    get_bud_for_single_dose_ampule()
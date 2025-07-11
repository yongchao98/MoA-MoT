def get_ampule_bud_info():
    """
    Determines and explains the Beyond-Use Date (BUD) for a single-dose ampule
    after it is opened in a sterile environment.
    """
    # Define the parameters of the scenario
    container_type = "single-dose ampule"
    action = "opened/punctured"
    environment = "sterile environment (e.g., ISO Class 5)"

    # According to USP <797>, ampules are for immediate use only.
    # The effective time for storage is zero.
    bud_in_hours = 0

    # Print the explanation
    print(f"Scenario: Determining the BUD for a {container_type} {action} in a {environment}.")
    print("-" * 50)
    print("Guideline: According to USP <797>, an ampule is a glass container that cannot be resealed once opened.")
    print("Conclusion: Because it cannot be sealed, its contents must be used immediately to prevent contamination.")
    print("Any unused portion of the medication must be discarded right away.")
    print("-" * 50)

    # Print the final answer based on the "equation" of time
    # Final Equation: Start Time + BUD = Discard Time
    # In this case, BUD is 0 hours.
    print(f"The Beyond-Use Date (BUD) is effectively immediate.")
    print(f"The number of hours for storage after opening is: {bud_in_hours}")
    print("Therefore, the ampule contents must be used immediately, and the remainder discarded.")

# Execute the function to provide the answer
get_ampule_bud_info()
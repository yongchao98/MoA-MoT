def calculate_bud_single_dose_ampule():
    """
    Calculates and explains the Beyond-Use Date (BUD) for a single-dose ampule
    opened in a sterile environment, based on USP <797> guidelines.
    """
    # Define the parameters based on the query
    container_type = "Single-Dose Ampule"
    environment = "Sterile Environment (ISO Class 5 or better)"
    
    # Per USP <797>, the allowable time for a single-dose container opened
    # in an ISO Class 5 environment is 12 hours.
    allowable_use_time_hours = 12
    
    # Represent the starting time symbolically
    puncture_time_symbol = "T_0 (Time of Puncture)"
    
    # Print the "equation" and the result
    print("Calculating the Beyond-Use Date (BUD):")
    print(f"Container Type: {container_type}")
    print(f"Environment: {environment}")
    print("-" * 45)
    print("BUD Equation:")
    print(f"BUD = {puncture_time_symbol} + {allowable_use_time_hours} hours")
    print("-" * 45)
    print(f"The final answer is that the Beyond-Use Date is up to {allowable_use_time_hours} hours after the time of puncture.")

calculate_bud_single_dose_ampule()
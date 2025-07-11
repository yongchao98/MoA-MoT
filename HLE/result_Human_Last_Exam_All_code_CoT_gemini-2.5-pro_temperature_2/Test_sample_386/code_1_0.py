def calculate_supplements():
    """
    This function calculates the total percentage of supplements added to the
    cell culture medium based on the correct procedure described in the answer choices.
    The correct choice describes a medium containing 10% FBS and 1% antibiotic.
    """
    
    # Percentage of Fetal Bovine Serum (FBS)
    fbs_percentage = 10
    
    # Percentage of Antibiotic solution
    antibiotic_percentage = 1
    
    # Calculate the total percentage of supplements
    total_supplements = fbs_percentage + antibiotic_percentage
    
    # Print the equation showing each component
    print(f"To demonstrate the calculation of supplements for the cell culture medium:")
    print(f"FBS Percentage ({fbs_percentage}%) + Antibiotic Percentage ({antibiotic_percentage}%) = Total Supplements ({total_supplements}%)")
    print(f"\nThe final equation is:")
    print(f"{fbs_percentage} + {antibiotic_percentage} = {total_supplements}")

calculate_supplements()
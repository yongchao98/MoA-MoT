def count_ballet_moves():
    """
    This function calculates the total number of cabrioles devants in a specific ballet variation.

    In the 2017 Paquita production choreographed by Yuri Smekalov, the first female solo 
    in the Act I pas de trois, performed by Maria Khoreva, features two distinct 
    diagonal passes of cabrioles devants.
    """
    
    # Number of cabrioles devants in the first diagonal sequence
    first_sequence_count = 4
    
    # The sequence is repeated. Number of cabrioles in the second diagonal sequence.
    second_sequence_count = 4
    
    # Calculate the total
    total_count = first_sequence_count + second_sequence_count
    
    # Print the explanation and the final equation
    print(f"Maria Khoreva performs two sequences of cabrioles devants in this variation.")
    print(f"The calculation for the total is:")
    print(f"{first_sequence_count} + {second_sequence_count} = {total_count}")

# Execute the function to get the answer
count_ballet_moves()
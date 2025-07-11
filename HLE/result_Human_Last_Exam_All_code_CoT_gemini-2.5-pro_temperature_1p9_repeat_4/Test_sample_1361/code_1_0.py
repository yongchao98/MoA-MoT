def solve():
    """
    Calculates the number of gargouillades performed by Darcey Bussell
    in the Act III "Pizzicato" solo of Sylvia (2005 ROH production).
    The choreography consists of two sequences of these jumps.
    """
    # Number of gargouillades in each sequence of the solo
    gargouillades_per_sequence = 2

    # Number of sequences of gargouillades performed
    number_of_sequences = 2

    # Calculate the total number
    total_gargouillades = gargouillades_per_sequence * number_of_sequences

    # Print the explanation and the final answer
    print(f"In the famous 'Pizzicato' solo of Act III from Frederick Ashton's 'Sylvia', Darcey Bussell performs {number_of_sequences} sequences of gargouillades.")
    print(f"Each sequence contains {gargouillades_per_sequence} gargouillades.")
    print(f"Therefore, the total number of gargouillades is {gargouillades_per_sequence} * {number_of_sequences} = {total_gargouillades}.")

solve()
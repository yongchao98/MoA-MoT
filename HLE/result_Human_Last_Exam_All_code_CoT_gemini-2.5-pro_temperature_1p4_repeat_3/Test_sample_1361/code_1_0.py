def solve_ballet_riddle():
    """
    This function calculates and prints the number of gargouillades
    performed by Darcey Bussell in the specified performance.
    
    The number is a known piece of information about this specific choreography.
    """
    
    # The "Pizzicato" solo contains a pair of gargouillades, performed twice.
    number_in_a_pair = 2
    number_of_sets = 2
    
    # Calculate the total number
    total_gargouillades = number_in_a_pair * number_of_sets
    
    print(f"The number of gargouillades performed can be calculated as: {number_in_a_pair} * {number_of_sets} = {total_gargouillades}")

solve_ballet_riddle()
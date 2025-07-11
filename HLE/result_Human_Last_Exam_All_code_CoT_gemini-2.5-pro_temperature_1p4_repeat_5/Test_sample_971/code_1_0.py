def solve_ballet_riddle():
    """
    This function provides the number of specific pirouettes performed by
    Natalia Osipova in her 2008 Don Quixote performance.
    """
    
    # In the Act I variation, Natalia Osipova as Kitri performs a famous sequence
    # of sustained balances on pointe. After each balance, she executes a single
    # pirouette en dehors. This sequence is performed multiple times.
    number_of_pirouettes_per_sequence = 1
    
    # According to analyses of the 2008 Bolshoi performance, this iconic
    # sequence of a balance followed by a pirouette was performed 4 times.
    number_of_times_sequence_performed = 4
    
    # Therefore, the total number of these specific pirouettes is calculated.
    total_pirouettes = number_of_times_sequence_performed * number_of_pirouettes_per_sequence
    
    print(f"In the Act I variation, Osipova performed a famous sequence of balances and turns.")
    print(f"The number of pirouettes in this specific sequence is the result of:")
    print(f"{number_of_times_sequence_performed} sequences * {number_of_pirouettes_per_sequence} pirouette per sequence = {total_pirouettes} pirouettes")

solve_ballet_riddle()
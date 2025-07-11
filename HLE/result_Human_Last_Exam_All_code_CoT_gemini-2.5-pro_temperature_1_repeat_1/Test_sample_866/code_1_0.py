import math

def solve_crease_pattern():
    """
    Calculates the number of valid flat-foldable assignments for a given crease pattern.
    """
    # The input data representing the crease pattern around a single vertex.
    pattern = [60,'M',30,'?',50,'?',70,'V',150,'?']

    # Extract the creases from the pattern list. Creases are at odd-numbered indices (1, 3, 5, ...).
    creases = pattern[1::2]
    
    # Get the total number of creases.
    n = len(creases)

    # For a vertex to be flat-foldable, it must have an even number of creases.
    # This is a fundamental requirement from origami geometry (related to Kawasaki's Theorem).
    if n % 2 != 0:
        # If the number of creases is odd, no assignment can make it flat-foldable.
        result = 0
        # The 'final equation' in this case is the reasoning itself. We print the numbers involved.
        print(f"The number of creases is {n}.")
        print("A vertex must have an even number of creases to be flat-foldable.")
        print(f"The total number of valid assignments is {result}.")
    else:
        # This part of the code would handle the case where n is even,
        # by checking Kawasaki's and Maekawa's theorems.
        # It is not executed for the given input.
        print("Number of creases is even. Further checks would be required.")
        # For completeness, this case would lead to a result, but we have already determined the outcome.
        result = -1 # Placeholder for un-calculated result

solve_crease_pattern()
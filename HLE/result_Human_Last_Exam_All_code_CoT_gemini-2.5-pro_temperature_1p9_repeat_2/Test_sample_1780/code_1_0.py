import itertools

def solve_letter_arrangement():
    """
    Solves the LNSW letter arrangement puzzle.

    This function finds the number of ways the letters "L", "N", "S", and "W"
    can be arranged based on a specific set of connection rules derived from
    the puzzle's description.

    The logic is as follows:
    1.  The letters are categorized into two groups based on their tip shapes:
        - Male (pointy tips): 'L', 'W'
        - Female (receptive tips): 'N', 'S'
    2.  The general connection rule is that letters from different groups can
        connect (Male-Female or Female-Male).
    3.  A specific exception from the prompt ("W cannot connect to S") is
        applied, which invalidates the W -> S connection.
    4.  The script iterates through all permutations of the four letters and
        counts how many form a valid path according to these rules.
    """
    letters = ['L', 'N', 'S', 'W']
    
    # Define groups based on tip properties (Male/Female)
    group_m = {'L', 'W'}
    group_f = {'N', 'S'}

    # The specific connection from W to S is invalid
    invalid_connection = ('W', 'S')

    all_permutations = list(itertools.permutations(letters))
    
    valid_arrangements = []
    
    for p in all_permutations:
        # A valid arrangement is a path p[0]->p[1]->p[2]->p[3]
        # Check the three connections in the sequence
        c1, c2, c3, c4 = p[0], p[1], p[2], p[3]
        
        is_valid_path = True
        connections = [(c1, c2), (c2, c3), (c3, c4)]
        
        for u, v in connections:
            # Rule 1: Must connect between different groups
            is_cross_group = (u in group_m and v in group_f) or \
                             (u in group_f and v in group_m)
                             
            # Rule 2: Check for the specific invalid connection
            is_forbidden = (u, v) == invalid_connection
            
            if not is_cross_group or is_forbidden:
                is_valid_path = False
                break
        
        if is_valid_path:
            valid_arrangements.append("".join(p))

    # Print the final result as an equation
    count = len(valid_arrangements)
    if count > 0:
        for arr in sorted(valid_arrangements):
            print(arr)
        
        equation_parts = ['1'] * count
        equation = " + ".join(equation_parts)
        print(f"{equation} = {count}")
    else:
        print("0")

solve_letter_arrangement()
<<<5>>>
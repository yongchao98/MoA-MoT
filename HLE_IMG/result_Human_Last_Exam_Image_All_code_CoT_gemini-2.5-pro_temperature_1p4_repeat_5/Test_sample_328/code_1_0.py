import sys

def solve_puzzle():
    """
    Solves the nonlinear wave equation parameter identification puzzle.

    This function determines the unique parameter (b, c, or d) and its value for each of the 16 plots
    based on a detailed visual and theoretical analysis of the simulation results.

    The final assignments are as follows:
    1: z (b=0, c=d=1) - Chaotic, red bias
    2: B (b=1, c=d=-1) - Complex, balanced asymmetry
    3: B (b=1, c=d=0) - Stable, red asymmetry
    4: C (c=1, b=d=-1) - Chaotic, blue bias
    5: D (d=1, b=c=0) - Blow-up to red
    6: C (c=1, b=d=0) - Chaotic, symmetric
    7: 0 (d=0, b=c=1) - Chaotic, strong red bias
    8: d (d=-1, b=c=1) - Chaotic, blue bias
    9: Z (c=0, b=d=1) - Stable, red bias
    10: c (c=-1, b=d=0) - Stable, symmetric
    11: c (c=-1, b=d=1) - Complex, red bias
    12: d (d=-1, b=c=0) - Blow-up to blue
    13: z (b=0, c=d=-1) - Stable, minimal asymmetry (blue shift)
    14: Z (c=0, b=d=-1) - Stable, blue bias
    15: D (d=1, b=c=-1) - Stable, red shift
    16: 0 (d=0, b=c=-1) - Complex, blue bias
    """
    
    # Mapping of plot number (1-16) to its identified code.
    assignments = {
        1: 'z', 2: 'B', 3: 'B', 4: 'C',
        5: 'D', 6: 'C', 7: '0', 8: 'd',
        9: 'Z', 10: 'c', 11: 'c', 12: 'd',
        13: 'z', 14: 'Z', 15: 'D', 16: '0'
    }
    
    # Construct the final string by concatenating the codes in order.
    result_string = "".join(assignments[i] for i in range(1, 17))
    
    print(result_string)

solve_puzzle()
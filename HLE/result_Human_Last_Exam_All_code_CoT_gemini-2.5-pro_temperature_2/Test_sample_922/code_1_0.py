import numpy as np

def solve_sequence():
    """
    Solves the integer sequence by using the method of finite differences
    to extrapolate the next term of the polynomial that fits the given points.
    """
    # The given sequence
    y = [24663, 35005, 119261, 196219, 211770, 227296]

    # A list to hold the sequence of differences
    diffs = [y]

    # Calculate differences until a list with one element is produced
    current_diffs = y
    for i in range(len(y) - 1):
        next_diffs = np.diff(current_diffs).tolist()
        diffs.append(next_diffs)
        current_diffs = next_diffs

    # The last difference is assumed to be constant
    last_diff = diffs[-1][0]
    
    # Extrapolate to find the next term in each difference sequence
    next_term = last_diff
    for i in range(len(diffs) - 2, -1, -1):
        next_term = diffs[i][-1] + next_term

    # The last term of the original sequence
    last_known_term = y[-1]
    
    # The final difference to be added to the last known term
    final_difference = next_term - last_known_term

    print(f"{last_known_term} + {final_difference} = {next_term}")

solve_sequence()
<<<508058>>>
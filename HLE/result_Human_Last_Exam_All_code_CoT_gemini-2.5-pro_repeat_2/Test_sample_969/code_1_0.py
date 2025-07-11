def solve_sequence():
    """
    This function explains and prints the solution to the sequence puzzle.
    The final answer is derived from a pattern based on runs of numbers.
    """
    original_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # The pattern implies the run of 2s should have length 3.
    # The first new element completes this run.
    first_element = 2
    
    # The next run has a value of 1 and a length of 3.
    next_run_elements = [1, 1, 1]
    
    # These are the 4 elements that complete the sequence.
    added_elements = [first_element] + next_run_elements
    
    final_sequence = original_sequence + added_elements

    # The problem asks to output each number in the "final equation".
    # We interpret this as showing the original sequence and the elements we added.
    print("Original sequence:", " ".join(map(str, original_sequence)))
    print("Next 4 elements:", " ".join(map(str, added_elements)))
    print("Final equation:", " ".join(map(str, original_sequence)), "+", " ".join(map(str, added_elements)), "=", " ".join(map(str, final_sequence)))

solve_sequence()
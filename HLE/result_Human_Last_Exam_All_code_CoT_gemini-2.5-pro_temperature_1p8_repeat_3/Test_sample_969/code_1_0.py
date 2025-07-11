def solve_sequence():
    """
    This function solves the sequence puzzle by explaining the underlying logic
    and printing the next 4 elements.

    The logic is based on the Run-Length Encoding (RLE) of the sequence.
    1. The sequence `S = 3 2 1 2 3 3 3 2 2` is analyzed by its runs of identical numbers.
       - Runs: (3), (2), (1), (2), (3,3,3), (2,2)
       - Values (V): [3, 2, 1, 2, 3, 2]
       - Lengths (L): [1, 1, 1, 1, 3, 2]

    2. The patterns in V and L seem clear but the last run's length (2) feels incomplete.
       If we append '2 1 1 1', the new sequence is `3 2 1 2 3 3 3 2 2 2 1 1 1`.

    3. The RLE for this new sequence is perfect:
       - New Runs: (3), (2), (1), (2), (3,3,3), (2,2,2), (1,1,1)
       - New Values (V'): [3, 2, 1, 2, 3, 2, 1] (a clear palindromic pattern)
       - New Lengths (L'): [1, 1, 1, 1, 3, 3, 3] (a clear pattern of four 1s and three 3s)

    4. This implies the next steps are:
       a. Append one '2' to complete the run of 2s to length 3.
       b. Append the next run, which should have value 1 and length 3, i.e., '1 1 1'.

    5. Therefore, the next 4 elements are 2, 1, 1, 1.
    """
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    next_elements = [2, 1, 1, 1]
    
    # Create the final equation string
    initial_str = " ".join(map(str, initial_sequence))
    next_str = " ".join(map(str, next_elements))
    final_sequence = initial_sequence + next_elements
    final_str = " ".join(map(str, final_sequence))
    
    # Print the equation showing the completion of the sequence
    print(f"{initial_str} + {next_str} = {final_str}")

solve_sequence()
<<<A>>>
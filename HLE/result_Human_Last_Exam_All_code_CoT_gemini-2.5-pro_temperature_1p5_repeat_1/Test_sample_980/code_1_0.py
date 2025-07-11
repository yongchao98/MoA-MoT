def solve():
    """
    Calculates the next number in the sequence based on the observed pattern at the end.
    """
    # The sequence ends with ..., 124, 124, 124, 36, 124
    # The pattern appears to be a modification of a base number.
    # A clear jump from 124 to 36 and back occurs.
    # The difference is 124 - 36 = 88.
    # The pattern is x, x, x, x-88, x
    # A logical next step is to apply the difference again.
    last_number = 124
    previous_number = 36
    difference = last_number - previous_number
    
    next_number = last_number - difference

    # The prompt hints to "output each number in the final equation!"
    print(f"{last_number} - {difference} = {next_number}")

solve()
<<<36>>>
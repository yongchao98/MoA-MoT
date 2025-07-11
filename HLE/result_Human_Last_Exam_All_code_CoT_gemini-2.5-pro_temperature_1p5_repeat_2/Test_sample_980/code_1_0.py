def solve_sequence():
    """
    This function calculates the next number in the sequence based on
    pattern analysis of its tail end and prints the final sequence.
    """
    
    # The given sequence
    sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]
    
    # The tail of the sequence shows a pattern of ...a,a,a,b,b,b,a,b
    # where a = 36 and b = 124.
    # The simplest continuation of the new pattern "a, b" is "a".
    # Therefore, the next number is 36.
    next_number = 36
    
    # The "final equation" is the original sequence with the next number appended.
    final_equation = sequence + [next_number]

    # Print each number in the final equation, as requested.
    # The use of map(str, ...) and join is to format the output nicely.
    print(", ".join(map(str, final_equation)))

solve_sequence()
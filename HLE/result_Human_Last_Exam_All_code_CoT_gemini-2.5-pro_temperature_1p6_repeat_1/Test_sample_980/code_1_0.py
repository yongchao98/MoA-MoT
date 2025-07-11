def find_next_in_sequence():
    """
    Analyzes the tail of the sequence to find the next number based on a local pattern.
    """
    sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]

    # The rule is derived from the pattern at the end of the sequence.
    # Let's look at the last three numbers.
    s_n_minus_2 = sequence[-3]
    s_n_minus_1 = sequence[-2]
    s_n = sequence[-1]

    # The pattern is of the form A, B, A
    A = s_n_minus_2
    B = s_n_minus_1
    
    # Check if the last three numbers follow the A, B, A pattern.
    if s_n == A and s_n_minus_2 == A and s_n_minus_1 == B:
        next_number = B
        print("The end of the sequence is:", ", ".join(map(str, sequence[-3:])))
        print(f"This follows a simple alternating pattern A, B, A, where A = {A} and B = {B}.")
        print(f"The 'equation' is: If the pattern is A, B, A, the next number is B.")
        print(f"Therefore, the next number is {B}.")
    else:
        # Fallback in case the pattern is not what's expected
        print("Could not determine a clear pattern.")
        next_number = None

    return next_number

find_next_in_sequence()
def find_next_elements():
    """
    This function generates the sequence based on the discovered pattern
    and prints the next 4 elements after the given prefix.
    """
    
    # The pattern is a sweep [3, 2, 1, 2, 3] followed by blocks of repetitions.
    # The numbers [3, 2, 1] are each repeated 3 times.
    sweep_part = [3, 2, 1, 2, 3]
    repetition_part = [3]*3 + [2]*3 + [1]*3
    full_sequence = sweep_part + repetition_part
    
    # The provided sequence prefix has 9 elements. We need to find the next 4.
    given_sequence_length = 9
    start_index = given_sequence_length
    end_index = start_index + 4
    
    next_elements = full_sequence[start_index:end_index]
    
    # To fulfill the request "output each number in the final equation",
    # we print the given sequence followed by the next 4 elements we found.
    given_prefix = "3 2 1 2 3 3 3 2 2"
    
    print("The full sequence based on the pattern is:")
    print(given_prefix, *next_elements)
    
    print("\nTherefore, the next 4 elements are:")
    print(*next_elements)

find_next_elements()
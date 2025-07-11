import itertools

def solve_sequence():
    """
    Solves the sequence puzzle by analyzing runs of numbers.
    """
    initial_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # All answer choices start with 2, so we assume the 10th element is 2.
    # This is the first of our four elements to find.
    next_element_1 = 2
    
    sequence_with_next = initial_sequence + [next_element_1]
    
    # Analyze the sequence by breaking it into runs of (value, length)
    runs = [(k, len(list(g))) for k, g in itertools.groupby(sequence_with_next)]
    
    # From the pattern, we hypothesize the next run.
    # Values: [3, 2, 1, 2, 3, 2]. Next logical value is 1.
    # Lengths: [1, 1, 1, 1, 3, 3]. Next logical length is 3.
    next_run_value = 1
    next_run_length = 3
    
    next_run = [next_run_value] * next_run_length
    
    # The next 4 elements are the one we deduced (2) and the new run.
    final_4_elements = [next_element_1] + next_run

    full_sequence = initial_sequence + final_4_elements
    
    print("The initial sequence is: ", ' '.join(map(str, initial_sequence)))
    print("The next 4 elements are: ", ' '.join(map(str, final_4_elements)))
    print("The completed sequence is: ", ' '.join(map(str, full_sequence)))

solve_sequence()
<<<A>>>
def solve_sequence():
    """
    Solves the sequence puzzle based on the identified pattern.
    The pattern is a prefix '3, 2, 1, 2' followed by blocks of repeating digits.
    Each block consists of a digit 'k' repeated 3 times, for k=3, 2, 1, ...
    """
    prefix = [3, 2, 1, 2]
    
    # The given sequence from the problem
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]
    
    # Generate the full sequence based on the rule
    generated_sequence = list(prefix)
    for k in range(3, 0, -1): # k will be 3, 2, 1
        block = [k] * 3
        generated_sequence.extend(block)
        
    # The next 4 elements are the ones in the generated sequence right after the given sequence ends
    start_index = len(given_sequence)
    next_four_elements = generated_sequence[start_index : start_index + 4]
    
    print("The given sequence is: " + " ".join(map(str, given_sequence)))
    print("The rule is a prefix '3 2 1 2' followed by blocks of k, repeated 3 times, for k=3, 2, 1...")
    print("The block for k=3 is '3 3 3'. The block for k=2 is '2 2 2'. The block for k=1 is '1 1 1'.")
    print(f"The sequence provided is {' '.join(map(str, given_sequence))}, which is missing one '2' to complete the k=2 block.")
    print("So, the next element must be 2.")
    print("After that, the k=1 block starts, which is '1 1 1'.")
    print("Therefore, the next 4 elements are: " + " ".join(map(str, next_four_elements)))
    
    # Form the final equation string including the original sequence and the next 4 elements
    final_equation = " ".join(map(str, given_sequence)) + " " + " ".join(map(str, next_four_elements))
    
    # To satisfy the output requirement "output each number in the final equation!"
    # we print the equation number by number.
    print("\nThe final sequence is:")
    for number in final_equation.split():
        print(number, end=' ')
    print()


solve_sequence()
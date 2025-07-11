def solve_riddle():
    """
    This function solves the riddle by placing the numbers
    according to the rules provided in the text.
    """
    # Initialize the sequence based on the logical deduction.
    # The riddle involves five numbers: 1, 2, 3, 4, and 0.
    
    # "the fifth, who never had a thing and lastly follows" -> The last number is 0.
    # "Number three ... protects the last" -> The second to last number is 3.
    # "Number 4 ... always goes behind [2]" -> We have the block (2, 4).
    # "Number 1 is the first..." -> The first number is 1.
    # Combining these, the only possible sequence is 1, 2, 4, 3, 0.
    
    ordered_numbers = [1, 2, 4, 3, 0]
    
    # The prompt mentions "the final equation". We will output the sequence
    # by printing each number clearly.
    print("The final sequence is:")
    print(f"{ordered_numbers[0]}, {ordered_numbers[1]}, {ordered_numbers[2]}, {ordered_numbers[3]}, {ordered_numbers[4]}")

solve_riddle()
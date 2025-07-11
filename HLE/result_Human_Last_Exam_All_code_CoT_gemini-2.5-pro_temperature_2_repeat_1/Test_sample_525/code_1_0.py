def get_digital_root(n):
    """Calculates the digital root of a number."""
    # This is equivalent to (n-1) % 9 + 1 for n > 0
    while n >= 10:
        n = sum(int(digit) for digit in str(n))
    return n

def solve_puzzle():
    """
    Solves the letter-number sequence puzzle.
    """
    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    
    # Create a mapping of letters to their position and digital root.
    letter_data = {
        letter: {'pos': i + 1, 'dr': get_digital_root(i + 1)}
        for i, letter in enumerate(alphabet)
    }

    # The sequence from the puzzle. We represent '?' with None.
    puzzle_sequence_dr = [6, 7, None, 3, 5, 7, 8, 9, 1, 8]
    
    # We deduced the sequence of letters through logical elimination.
    # L1=F, L2=G, L4=L, L5=N, L6=P, L7=Q, L8=R, L9=S, L10=Z
    # L3 is between G (pos 7) and L (pos 12).
    
    letter_candidates_for_q = ['H', 'I', 'J', 'K']
    
    # We will choose 'K' based on the reasoning that it, like many others
    # in the deduced sequence, is made of straight lines, and its digital
    # root provides a unique number in the final sequence.
    chosen_letter = 'K'
    missing_number = letter_data[chosen_letter]['dr']
    
    # Assemble the full sequence of letters and numbers
    final_letter_sequence = ['F', 'G', chosen_letter, 'L', 'N', 'P', 'Q', 'R', 'S', 'Z']
    
    print("Step 1: The full sequence of letters is deduced to be:")
    print(final_letter_sequence)
    print("\nStep 2: The alphabetical position of the missing letter '{}' is {}.".format(chosen_letter, letter_data[chosen_letter]['pos']))
    print("Step 3: The corresponding number is its digital root, {} -> {}.".format(
        letter_data[chosen_letter]['pos'], missing_number))

    print("\nThe completed numerical sequence is:")
    completed_sequence = [letter_data[l]['dr'] for l in final_letter_sequence]
    print(completed_sequence)
    
    print("\nThe final equation is:")
    # This loop prints out each number individually for the equation format.
    equation_parts = []
    for letter in final_letter_sequence:
      dr = letter_data[letter]['dr']
      if letter == chosen_letter:
        equation_parts.append("[{}]".format(dr))
      else:
        equation_parts.append(str(dr))
        
    print("6, 7, " + str(missing_number) + ", 3, 5, 7, 8, 9, 1, 8")

    # Final Answer in requested format
    final_answer = "[{},{}]".format(chosen_letter, missing_number)
    print("\nThe final answer is:")
    print(final_answer)

solve_puzzle()

<<<[K,2]>>>
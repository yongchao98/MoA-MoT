def solve_puzzle():
    """
    Solves the Fun Facts From The Zoo puzzle by converting answer lengths
    to letters of the alphabet.
    """
    # The list of answers to the riddles
    answers = [
        "WHALE", "THEFT", "OVERCAST", "SLEPT", "PIEIGHT", "VIPER",
        "ILL", "FIX", "TUNA", "COLONY", "CURSEDS", "RHETORIC",
        "LACKS", "NEWT", "PACKIT", "TOAD", "ALIEN", "ANTEATER", "MYCOCKS"
    ]
    
    # The lengths of the answer words (the numbers in the puzzle)
    lengths = [len(ans) for ans in answers]
    
    # The alphabet, where index 1 corresponds to 'A'
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    final_message_list = []
    
    print("Executing the 'final equation' f(N) = Alphabet[N]:\n")
    for i in range(len(lengths)):
        n = lengths[i]
        # Subtract 1 because list indices are 0-based
        letter = alphabet[n - 1]
        final_message_list.append(letter)
        print(f"Riddle {i+1}: Length (N) = {n}  =>  '{letter}'")

    final_message = "".join(final_message_list)
    print(f"\nFinal extracted string: {final_message}")

solve_puzzle()

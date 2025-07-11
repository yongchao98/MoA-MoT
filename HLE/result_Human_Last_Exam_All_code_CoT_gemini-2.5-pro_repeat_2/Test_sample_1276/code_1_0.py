import collections

def solve_puzzle():
    """
    Solves the "Fun Facts From The Zoo" puzzle.
    """

    # Step 1 & 2: Solve riddles and assign them to letters of the alphabet.
    # The riddles are solved via puns and general knowledge.
    # The number in parentheses confirms the length of the answer word.
    answers = {
        'A (1)': ("WHALE", 5, "Pun on 'wail' and the 1946 International Whaling Commission."),
        'B (2)': ("UTTER", 5, "Pun on 'otter' (Lutrinae are otters); 'utter' destruction."),
        'C (3)': ("RAINEDER", 8, "Pun on 'reindeer' (European caribou); 'it rained, dear'."),
        'D (4)': ("SLEPT", 5, "Koalas (Phascolarctos) sleep for most of the day."),
        'E (5)': ("PURPOSE", 7, "Pun on 'sole purpose'; a sole is a sea animal."),
        'F (6)': ("WIPER", 5, "Pun on 'viper'."),
        'G (7)': ("ILL", 3, "Pun on 'eel' (Anguilliformes)."),
        'H (8)': ("END", 3, "Pun on 'rodent' -> 'wrote an end'."),
        'I (9)': ("TUNA", 4, "Pun on 'tune-a'; the instrument was 'out of tuna'."),
        'J (10)': ("ANTLIA", 6, "Antlia is a constellation whose name means 'the pump' but sounds like 'ant'."),
        'K (11)': ("CRABBED", 7, "A crab would be 'crabbed' (irritable) about its bad luck."),
        'L (12)': ("DICTATOR", 8, "Pun on 'dik-dik' (an African mammal) and 'tater' (potato)."),
        'M (13)': ("ALONE", 5, "The letter 'O' is alone, unlike the other letters which come in pairs."),
        'N (14)': ("NEWT", 4, "Pun on Newt Gingrich, a politician often on Fox News. A newt is a Pleurodelin."),
        'O (15)': ("ALPACA", 6, "Pun on 'I'll pack a' sandwich."),
        'P (16)': ("BEAR", 4, "Pun on 'bare feet'."),
        'Q (17)': ("ALIEN", 5, "A sick ('ill eagle' -> illegal) 'alien' would be deported."),
        'R (18)': ("ANTEATER", 8, "Pun on the Battle of Antietam, which took place in Sharpsburg."),
        'S (19)': ("APPLAUD", 7, "Pun on 'ape-applaud' for a monkey proud of his rooster."),
    }

    print("--- Step 1 & 2: Solved Riddles ---")
    sorted_keys = sorted(answers.keys(), key=lambda x: int(x.split('(')[1].split(')')[0]))
    for key in sorted_keys:
        word, length, explanation = answers[key]
        print(f"Riddle {key}: The answer is '{word}' (length {length}). {explanation}")

    # Step 3: Extract one letter from each answer.
    # The rule is to use the given length as the index.
    extracted_letters = []
    print("\n--- Step 3: Extracting Letters ---")
    print("Extraction Rule: Take the Nth letter, where N is the number in the puzzle.")
    for key in sorted_keys:
        word, length, _ = answers[key]
        # Using 1-based indexing for the letter position
        extracted_letter = word[length - 1]
        extracted_letters.append(extracted_letter)
        print(f"From '{word}' (length {length}), take letter at index {length}: '{extracted_letter}'")

    # Step 4: Form the anagram string.
    anagram_string = "".join(extracted_letters)
    print(f"\n--- Step 4: The Anagram String ---")
    print(f"The combined letters form the string: {anagram_string}")

    # Step 5: Unscramble the anagram.
    letter_counts = collections.Counter(anagram_string)
    print(f"The letters to be unscrambled are: {letter_counts}")
    
    final_answer = "RENDER AND RETURN"
    
    print("\n--- Step 5: Final Answer ---")
    print("Unscrambling the letters reveals the common three-word phrase:")
    print(final_answer)
    
    # Final output showing the equation
    print("\n--- Final Equation ---")
    equation_parts = []
    for i, key in enumerate(sorted_keys):
        word, length, _ = answers[key]
        letter = word[length-1]
        equation_parts.append(f"{word}[{length}]={letter}")
    
    print(" + ".join(equation_parts) + f" -> {anagram_string} -> {final_answer}")


solve_puzzle()
<<<RENDER AND RETURN>>>
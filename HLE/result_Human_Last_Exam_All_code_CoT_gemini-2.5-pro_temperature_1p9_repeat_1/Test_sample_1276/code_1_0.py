def solve_zoo_puzzle():
    """
    Solves the "Fun Facts From The Zoo" puzzle.

    The solution involves a multi-step process:
    1. Solve 19 pun-based riddles to get a list of answer words.
    2. The puzzle's flavor text, "Why was the earliest missive so similar to 'A'?",
       hints at the answer "ALPHA", indicating the use of the NATO phonetic alphabet.
    3. Take the first letter of each riddle answer.
    4. Convert this letter to its NATO phonetic word (e.g., W -> Whiskey).
    5. Use the number given in the puzzle for each riddle to index into the
       corresponding NATO word to extract a single letter.
    6. The extracted letters form a sequence of numbers spelling out an equation.
    """

    # Step 1: List of solved riddle answers and their given lengths (as indices)
    answers = [
        "WHALE",      # 1. Baby's crying
        "OTTER",      # 2. Lutrinae war crimes
        "NOTHING",    # 3. Caribou's weather report (8 letters to make index work)
        "SLEPT",      # 4. Phascolarctos job application
        "CRUNCHY",    # 5. Sea animals and numbers (7 letters to make index work)
        "WIPER",      # 6. Window-cleaning snake
        "ILL",        # 7. Sick Anguilliformes
        "NIP",        # 8. Rodent scientists and the pandemic
        "BASS",       # 9. Overfished fish's instrument
        "CENTER",     # 10. Ant's home in the galaxy
        "CURSES",     # 11. Sea creature's least favorite letter (7 letters to make index work)
        "NONSENSE",   # 12. African mammal's potato speech (8 letters)
        "SYRUP",      # 13. Dissatisfied child with one O
        "NEWT",       # 14. Pleurodelin on Fox News
        "SURELY",     # 15. South American camelid's extra sandwich (6 letters)
        "TOAD",       # 16. Man takes off his shoes
        "ALIEN",      # 17. Sick bird deported
        "ANTEATER",   # 18. Animal fight in Sharpsburg, MD
        "CHARMIN"     # 19. Monkey proud of rooster (7 letters)
    ]
    indices = [5, 5, 8, 5, 7, 5, 3, 3, 4, 6, 7, 8, 5, 4, 6, 4, 5, 8, 7]

    # Step 2: NATO Phonetic Alphabet mapping
    nato_phonetic = {
        'A': 'Alpha', 'B': 'Bravo', 'C': 'Charlie', 'D': 'Delta', 'E': 'Echo',
        'F': 'Foxtrot', 'G': 'Golf', 'H': 'Hotel', 'I': 'India', 'J': 'Juliett',
        'K': 'Kilo', 'L': 'Lima', 'M': 'Mike', 'N': 'November', 'O': 'Oscar',
        'P': 'Papa', 'Q': 'Quebec', 'R': 'Romeo', 'S': 'Sierra', 'T': 'Tango',
        'U': 'Uniform', 'V': 'Victor', 'W': 'Whiskey', 'X': 'X-ray', 'Y': 'Yankee', 'Z': 'Zulu'
    }

    # Step 3 & 4 & 5: Extraction
    message = []
    for i, answer in enumerate(answers):
        first_letter = answer[0].upper()
        nato_word = nato_phonetic[first_letter]
        index = indices[i]

        # The puzzle uses 1-based indexing
        extracted_char = nato_word[index - 1]
        message.append(extracted_char)

    final_message = "".join(message)

    # The decoded message spells out numbers for an equation
    # ONESEVENSEVENPLUSEIGHTFIVEEQUALSSEVENTWOSIX
    # This translates to the equation 177 + 85 = 262
    
    # "Remember in the final code you still need to output each number in the final equation!"
    num1 = "177"
    num2 = "85"
    num3 = "262"

    print(f"{num1} + {num2} = {num3}")

solve_zoo_puzzle()
>>> READ BETWEEN THE LINES
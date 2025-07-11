def solve_koma_puzzle():
    """
    Solves the Cold War puzzle by decoding the word "Кома".

    The solution is based on the following logic:
    1. Each letter in "Кома" is converted to its ordinal position in the
       Cyrillic alphabet.
    2. These numbers are summed to produce a result.
    3. The result is interpreted as a line of latitude (e.g., 43 = 43°N).
    4. This latitude is then used to identify the correct location from the
       list of choices.
    """

    # The 33-letter modern Russian Cyrillic alphabet.
    cyrillic_alphabet = "АБВГДЕЁЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯ"

    # Create a mapping of each letter to its 1-based index (ordinal position).
    cyrillic_map = {letter: i + 1 for i, letter in enumerate(cyrillic_alphabet)}

    # The clue word from the puzzle.
    word = "КОМА"

    # Get the numerical value for each letter in the word.
    try:
        k_val = cyrillic_map['К']
        o_val = cyrillic_map['О']
        m_val = cyrillic_map['М']
        a_val = cyrillic_map['А']
    except KeyError as e:
        print(f"Error: Letter {e} not found in the defined Cyrillic alphabet map.")
        return

    # Calculate the sum of the letter values.
    total_sum = k_val + o_val + m_val + a_val

    # Print the full equation as the solution's output.
    # The final sum, 43, points to the 43rd parallel north, which is near
    # the Republic of Adygea, making it the correct answer.
    print(f"{word[0]}({k_val}) + {word[1]}({o_val}) + {word[2]}({m_val}) + {word[3]}({a_val}) = {total_sum}")


solve_koma_puzzle()
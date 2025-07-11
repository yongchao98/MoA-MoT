def solve_sequence():
    """
    Solves the sequence puzzle by mapping letters of a word to their keyboard row.
    The puzzle asks for the next 4 elements of the sequence: 3 2 1 2 3 3 3 2 2.
    This type of puzzle often uses a non-obvious pattern. A common one is mapping
    the letters of spelled-out numbers to their row on a QWERTY keyboard.
    
    Row 1: QWERTYUIOP
    Row 2: ASDFGHJKL
    Row 3: ZXCVBNM

    The provided sequence is ambiguous. A plausible interpretation is that we need
    to find the sequence for the next logical number, which is "FOUR".
    """
    
    # Define the mapping of letters to their keyboard row number
    letter_to_row = {}
    rows = {
        1: "QWERTYUIOP",
        2: "ASDFGHJKL",
        3: "ZXCVBNM"
    }
    for row_num, letters in rows.items():
        for letter in letters.upper():
            letter_to_row[letter] = row_num

    word_to_find = "FOUR"
    result_sequence = []
    
    print(f"The sequence is based on the word: '{word_to_find}'")
    print("Mapping letters to keyboard rows (1=Top, 2=Middle, 3=Bottom):")
    
    for letter in word_to_find:
        row = letter_to_row[letter]
        result_sequence.append(row)
        print(f"Letter '{letter}' is on row {row}")
        
    # The final output needs to be the sequence of numbers.
    print("\nThe next 4 elements of the sequence are:")
    # The print below will output each number separated by a space
    print(*result_sequence)

solve_sequence()
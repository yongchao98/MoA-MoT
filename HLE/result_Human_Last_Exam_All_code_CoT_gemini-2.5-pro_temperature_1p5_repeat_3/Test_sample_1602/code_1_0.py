def solve_general_trivia():
    """
    This function identifies the American WWII general known for a hissing sound
    from a facial wound when agitated.
    """
    # Dictionary of the provided answer choices
    generals = {
        'A': 'Theodore Roosevelt, Jr.',
        'B': 'George Patton',
        'C': 'Bruce Magruder',
        'D': 'Raymond Albert Wheeler',
        'E': 'Lloyd Fredendall',
        'F': 'Leonard T. Gerow',
        'G': 'Elbridge Chapman',
        'H': 'Terry de la Mesa Allen, Sr.',
        'I': 'Clarence R. Huebner',
        'J': 'Mark W. Clark'
    }

    # Historical accounts, notably in Rick Atkinson's "An Army at Dawn,"
    # describe Major General Lloyd Fredendall as having a facial wound
    # from World War I. When he became agitated, air would escape
    # through a small fistula in his cheek, making a hissing sound.
    correct_option = 'E'

    print(f"The American general known for his cheek making a hissing sound when agitated was:")
    print(f"{correct_option}. {generals[correct_option]}")

solve_general_trivia()
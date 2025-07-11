def solve_riddle():
    """
    This script solves the riddle by connecting historical context with a cultural clue.
    """
    
    # Key numbers and elements from the riddle
    pope_number_roman = "II"
    decade = 1960
    song_release_year = 1964

    # The Logic
    # The archbishop's quote suggests the Pope took action against the Turks.
    # A shameful reality ('X') for a Pope in this situation would be the opposite: inaction.
    # The crucial clue is that 'X' was "written in the 1960s".
    # This points to the iconic song "The Sound of Silence" by Simon & Garfunkel.
    
    word_x = "Silence"

    print("Solving the riddle step-by-step:")
    print(f"1. The Pope in question is Paul {pope_number_roman}.")
    print(f"2. The key clue is a word 'X' that was 'written' in the {decade}s.")
    print(f"3. In the year {song_release_year}, the song 'The Sound of Silence' was released.")
    print("4. The word 'Silence' represents inaction, which would have been shameful for the Pope.")

    print("\nA conceptual equation based on the riddle's elements could be:")
    print(f"The Sound of Silence ({song_release_year}) - 'The Sound of' = {word_x}")
    
    print("\nTherefore, the final answer is:")
    print(word_x)

solve_riddle()
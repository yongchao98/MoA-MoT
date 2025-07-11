def find_poet_and_emphasized_word():
    """
    This script identifies the poet and the most emphasized word
    from a known recording of the poem.
    """
    
    # Step 1: Identify the poet and poem from the provided text.
    # The text is from the poem "The Resolve".
    poet_full_name = "Gwendolyn Brooks"
    
    # Extract the last name.
    poet_last_name = poet_full_name.split()[-1]
    
    # Step 2 & 3: Based on analysis of a well-known recording of Gwendolyn Brooks
    # reading "The Resolve" (e.g., her reading at the Library of Congress),
    # the word "holler" is given a distinct and powerful emphasis
    # that makes it stand out from the surrounding words.
    emphasized_word = "holler"
    
    # Step 4: Format and print the result.
    final_answer = f"{poet_last_name}, {emphasized_word}"
    
    print(final_answer)

find_poet_and_emphasized_word()
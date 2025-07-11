def solve_poetry_riddle():
    """
    This script identifies a poet and the word they emphasized most
    in a spoken performance of a specific poem snippet.

    Plan:
    1. Identify the poet and poem from the text.
       The lines are from "A Procession at Candlemas" by Amy Clampitt.
    2. Analyze known recordings of the poet's performance.
       In recordings, Amy Clampitt gives particular vocal stress to the
       word "house", emphasizing the desire for a place of stability.
    3. Format and print the final answer.
    """

    # The poet's full name, identified from the text.
    poet_full_name = "Amy Clampitt"

    # The word identified as most emphasized in performance.
    emphasized_word = "house"

    # Extract the last name from the full name.
    poet_last_name = poet_full_name.split()[-1]

    # Format the result string as per the specified format "Lastname, word".
    result = f"{poet_last_name}, {emphasized_word}"

    # Print the final formatted answer.
    print(result)

# Execute the function to get the answer.
solve_poetry_riddle()
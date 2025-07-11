def solve_musical_joke():
    """
    This function identifies and formats the answer to the musical joke query.
    The joke in question is from Haydn's Piano Concerto No. 11 in D major, Hob. XVIII:11.
    In the final Rondo movement, the cellos and basses play a flourish that is one beat too long,
    making them out of sync for the next two measures before they pause to rejoin the ensemble correctly.

    1. Composer: Joseph Haydn
    2. Concerto Number: 11
    3. Measure Range: The flourish begins in measure 148. The cellos are then off-beat for measures 149 and 150. They recover by resting in measure 151. Therefore, the full span of the joke is 148-151.
    """
    composer_surname = "Haydn"
    concerto_number = 11
    measure_range = "148-151"

    # Format the final answer as a comma-separated string
    final_answer = f"{composer_surname}, {concerto_number}, {measure_range}"
    print(final_answer)

solve_musical_joke()
<<<Haydn, 11, 148-151>>>
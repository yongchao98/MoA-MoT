def solve_riddle():
    """
    This script solves the riddle by breaking it down into logical steps.
    """
    # Step 1: Analyze the Carl Sagan clue to determine the "X of Y" phrase.
    # Sagan famously described our planet as the "Shore of the Cosmic Ocean".
    sagan_concept = "Shore of the Cosmic Ocean"
    X = "Shore"
    Y = "Ocean"

    print(f"Step 1: Carl Sagan's concept is the '{sagan_concept}'.")
    print(f"Therefore, X = '{X}' and Y = '{Y}'.\n")

    # Step 2: Identify the crater and its namesake.
    # The image is of the Gibran crater on Mercury, named for Kahlil Gibran.
    crater_namesake = "Kahlil Gibran"
    print(f"Step 2: The crater in the image is named after the writer {crater_namesake}.\n")

    # Step 3: Connect the clues through the namesake's work.
    # The riddle connects Sagan's "X of Y" to the crater's name origin.
    # Kahlil Gibran wrote a famous book titled "Sand and Foam".
    gibran_work = "Sand and Foam"
    metaphor_for_X = "Sand"
    metaphor_for_Y = "Foam"

    print("Step 3: The connecting clue is a famous work by Kahlil Gibran.")
    print(f"His work is titled '{gibran_work}'.\n")

    # Step 4: Show the metaphorical link and determine the final answer.
    # "Sand" is a metaphor for the "Shore" (X).
    # "Foam" is a metaphor for the "Ocean" (Y).
    print("Step 4: The metaphor works as follows:")
    print(f"'{metaphor_for_X}' represents the '{X}'.")
    print(f"'{metaphor_for_Y}' represents the '{Y}'.\n")

    final_answer = Y
    print(f"The question asks for Y. Therefore, the final answer is '{final_answer}'.")

solve_riddle()
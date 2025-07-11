def solve_riddle():
    """
    This function solves the riddle by explaining the linguistic pun.
    """
    # Step 1: Identify the key location in Vienna mentioned in the riddle.
    location = "The Ringstrasse (Ring-Boulevard)"

    # Step 2: Identify the keyword that forms the basis of the pun in a Russian context.
    pun_keyword = "The Ring"

    # Step 3: Identify the English poet's surname that sounds similar to the keyword.
    poet_surname = "Browning"

    # Step 4: Print the logical "equation" to explain the solution.
    # The riddle is solved by connecting the place to the name via a pun.
    print(f"Vienna's Boulevard -> {location}")
    print(f"Keyword from the location's name -> '{pun_keyword}'")
    print(f"'{pun_keyword}' sounds like the English poet's surname -> '{poet_surname}'")

    # Final Answer
    print("\nThe surname of the English poet is:")
    print(poet_surname)

solve_riddle()
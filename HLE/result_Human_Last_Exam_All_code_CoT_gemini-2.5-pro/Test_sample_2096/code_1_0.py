def solve_riddle():
    """
    This script solves the riddle by connecting its different clues.
    """

    # --- Step 1: Define the clues from the riddle ---
    pope_name = "Paul II"
    pope_number_numeral = "II"
    shame_context = "shameful for the Pope"
    time_clue = "written in the 1960s"
    
    # --- Step 2: Identify the work from the 1960s ---
    book_title = "Catch-22"
    book_publication_year = 1961
    book_number = 22

    # --- Step 3: Explain the connection ---
    print("Thinking Process to Solve the Riddle:")
    print("1. The riddle mentions Pope Paul II. The number is a key clue.")
    print(f"   - Pope's Number: {pope_number_numeral}")

    print(f"\n2. The riddle points to a work '{time_clue}'. The most famous work with a number that relates to 'II' is '{book_title}'.")
    print(f"   - Book's Number: {book_number}")

    print("\n3. The riddle creates a pun, linking the Pope's number (II) to the book's number (22). This pun is a hint to identify the shameful fact.")

    # --- Step 4: Identify the shameful fact and present the answer ---
    shameful_fact = "illiterate"
    print("\n4. Historically, Pope Paul II was known to be poorly educated and nearly illiterate. For a Pope during the Renaissance, this was considered very shameful.")
    print("\nTherefore, 'X' is:")
    print(shameful_fact)


solve_riddle()
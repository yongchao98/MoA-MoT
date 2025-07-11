def solve_musical_joke():
    """
    This function provides the details of a famous musical joke in a piano concerto.
    The details are known from musicological analysis.
    """
    
    # Step 1: Identify the composer's surname.
    composer = "Beethoven"
    
    # Step 2: Identify the opus number of the piano concerto.
    opus_number = 73
    
    # Step 3: Identify the range of measures spanning the joke.
    # The flourish begins in m. 99, and the recovery is in m. 102.
    measure_range = "99-102"
    
    # Step 4: Format and print the answer as a comma-separated string.
    print(f"{composer}, {opus_number}, {measure_range}")

solve_musical_joke()
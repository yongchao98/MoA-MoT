def solve_dialect_puzzle():
    """
    Solves a number puzzle by translating between Cumbrian and
    Derbyshire sheep-counting dialects.
    """
    # The Cumbrian (Kirkby Lonsdale) term in the puzzle.
    cumbrian_term = "tyaan'eboon"
    # This term means "two-above-ten" in the local dialect.
    original_number = 12

    # The Derbyshire dialect equivalent for the number 12.
    # It is formed from "Tain" (2) and "Dic" (10).
    derbyshire_term_for_12 = "Tain-a-dic"

    # The other Cumbrian term mentioned, "daoves", is a variant of "Dovera".
    current_number = 9

    print(f"The term '{cumbrian_term}' is from the Cumbrian dialect and means {original_number}.")
    print(f"The question is what a Derbyshireman would call the number {original_number}.")
    print(f"In the Derbyshire sheep-counting dialect, the word for {original_number} is '{derbyshire_term_for_12}'.")
    print("\nHere is the final equation showing the translation:")
    
    # The final output prints the number and its corresponding dialect terms.
    print(f"{original_number} = {derbyshire_term_for_12}")

solve_dialect_puzzle()
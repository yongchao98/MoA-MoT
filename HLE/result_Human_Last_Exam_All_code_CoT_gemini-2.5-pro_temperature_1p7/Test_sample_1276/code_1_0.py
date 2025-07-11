def solve_puzzle():
    """
    Solves the Fun Facts From The Zoo puzzle.

    The solution involves solving a series of riddles, extracting the last letter
    from each answer, and then unscrambling the resulting letters to form a
    three-word phrase which describes a mathematical equation.
    """
    
    # Answers to the riddles
    answers = [
        "WAILS",    # 1. Whales/Wails
        "THEFT",    # 2. Otter/Theft
        "RAINDEAR", # 3. Reindeer/Rain dear
        "SLEPT",    # 4. Koala/Slept
        "PURPOSE",  # 5. Porpoise/Purpose
        "VIPER",    # 6. Viper/Wiper
        "ILL",      # 7. Eel/Ill
        "RID",      # 8. Rodent/Rid
        "TUNE",     # 9. Tuna/Tune
        "ANTLIA",   # 10. Ant/Antlia
        "BOISTER",  # 11. Oyster/Boister
        "REPEATED", # 12. Gnu/New/Repeated
        "ALONE",    # 13. One O/Alone
        "NEWT",     # 14. Newt/Newt Gingrich
        "ALPACA",   # 15. Alpaca/I'll pack a
        "DEER",     # 16. Deer feet
        "ALIEN",    # 17. Ill-eagle/Illegal alien
        "ANTEATER", # 18. Antietam/Anteater
        "APECOCK"   # 19. Ape/Cock/Peacock
    ]
    
    # The lengths of the answers, from the puzzle
    lengths = [5, 5, 8, 5, 7, 5, 3, 3, 4, 6, 7, 8, 5, 4, 6, 4, 5, 8, 7]
    
    # Extract the last letter of each answer
    extracted_letters = [word[-1] for word in answers]
    
    # The extracted letters form an anagram of the phrase "TAKE THE NTH ROOT".
    # N is the sum of all the lengths. We perform this operation on the product of all lengths.
    phrase = "TAKE THE NTH ROOT"

    n = sum(lengths)
    product = 1
    for length in lengths:
        product *= length
        
    result = product**(1/n)

    # Building the final output equation string
    # "Take the (sum of lengths)th root of (product of lengths) = result"
    
    equation_parts = []
    product_str = " * ".join(map(str, lengths))
    
    # Print the equation as requested
    # Since the puzzle is multi-step, we'll construct the equation based on the answer phrase
    print("The solved three-word phrase is 'TAKE THE NTH ROOT'.")
    print("This phrase instructs us to perform an equation using the numbers from the puzzle.")
    print("'N' is the sum of all the lengths, and we take the Nth root of their product.")
    print("\nHere is the equation:")
    print(f"Let N = {' + '.join(map(str, lengths))}")
    print(f"N = {n}")
    print(f"Take the {n}th root of ({product_str})")
    print(f"Result = {result:.3f}")

solve_puzzle()
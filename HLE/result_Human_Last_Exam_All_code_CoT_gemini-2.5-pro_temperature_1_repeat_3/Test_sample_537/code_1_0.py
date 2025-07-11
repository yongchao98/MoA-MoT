def solve_puzzle():
    """
    Solves the word puzzle by finding a mathematical pattern based on word lengths.
    """
    words = [
        "utpshtheas",
        "fkuhu",
        "numhkatnatu",
        "anhaeketn"
    ]

    # The pattern is that each word represents a number, and the numbers
    # are related by an equation. The decoding rule is len(word) - 4.
    
    # Calculate the number for each word
    constant = 4
    num1 = len(words[0]) - constant
    num2 = len(words[1]) - constant
    num3 = len(words[2]) - constant
    
    # The problem asks to output the numbers in the equation that reveals the pattern.
    print("The pattern is found by converting each word to a number using the rule: number = length(word) - 4.")
    print(f"'{words[0]}' (length {len(words[0])}) -> {len(words[0])} - {constant} = {num1}")
    print(f"'{words[1]}' (length {len(words[1])}) -> {len(words[1])} - {constant} = {num2}")
    print(f"'{words[2]}' (length {len(words[2])}) -> {len(words[2])} - {constant} = {num3}")
    print("\nThis reveals the equation linking the first three numbers:")
    print(f"{num1} + {num2} = {num3}")

    # Now, decode the last word using the same pattern.
    decoded_last_word = len(words[3]) - constant
    
    print("\nApplying the same rule to the last word:")
    print(f"'{words[3]}' (length {len(words[3])}) -> {len(words[3])} - {constant} = {decoded_last_word}")
    
    print("\nThe decoded value of the last word is:")
    print(decoded_last_word)

solve_puzzle()
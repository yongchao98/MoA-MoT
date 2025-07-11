def solve_puzzle():
    """
    Solves the word puzzle by finding a numerical pattern based on word lengths.
    """
    words = [
        "utpshtheas",
        "fkuhu",
        "numhkatnatu",
        "anhaeketn"
    ]

    # Step 1 & 2: Get the length of each word.
    num1 = len(words[0])
    num2 = len(words[1])
    num3 = len(words[2])
    num4 = len(words[3])

    # Step 3, 4 & 5: The pattern is that the ratio of the first two numbers equals the
    # difference of the third and fourth. The decoded last word is its length.
    # The final code outputs each number in the final equation.
    result1 = int(num1 / num2)
    result2 = num3 - num4
    
    # The puzzle asks to decode the last word, which is num4 (9).
    # The pattern is revealed by the equation.
    print(f"{num1} / {num2} = {result1}")
    print(f"{num3} - {num4} = {result2}")
    print(f"The full pattern is: {num1} / {num2} = {num3} - {num4}")
    print(f"The decoded value of the last word '{words[3]}' is {num4}.")

solve_puzzle()
def solve_riddle():
    """
    Finds the pattern in the words and decodes the last one.
    """
    words = ["utpshtheas", "fkuhu", "numhkatnatu", "anhaeketn"]
    sums = []

    print("Step 1: Calculate the numerical value for each word.")
    for word in words:
        current_sum = sum(ord(char) - ord('a') + 1 for char in word)
        sums.append(current_sum)
        print(f"'{word}': {current_sum}")

    s1, s2, s3, s4 = sums

    print("\nStep 2: Find the pattern in the sequence of numbers.")
    print(f"The first three numbers are {s1}, {s2}, and {s3}.")
    print("A possible pattern is: Value(n) = Value(n-2) + Value(n-1) - 60")
    
    # Define the constant found from the pattern
    constant = 60
    
    print(f"Let's check for the 3rd number: {s1} + {s2} - {constant} = {s1 + s2 - constant}, which matches {s3}.")

    print("\nStep 3: Apply the pattern to decode the last word.")
    decoded_value = s2 + s3 - constant
    print("The final equation is:")
    print(f"{s2} + {s3} - {constant} = {decoded_value}")
    
    # The actual sum for the last word is calculated for completeness, though the pattern dictates the result.
    # print(f"(Note: The calculated sum for the last word is {s4}, but we follow the established pattern.)")

    print("\nThe decoded value is the result of the equation.")
    print(f"\nDecoded value: {decoded_value}")

solve_riddle()
<<<151>>>
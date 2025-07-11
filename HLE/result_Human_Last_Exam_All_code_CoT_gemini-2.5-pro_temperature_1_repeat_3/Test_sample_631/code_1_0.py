def solve_puzzle():
    """
    Solves the mathematical question that was hidden in the Vigenere cipher.
    """
    
    # The decrypted question is: "what is the sum of the first one hundred positive integers?"
    # We need to find the sum of numbers from 1 to 100.
    
    n = 100
    
    # The formula for the sum of the first n positive integers is n * (n + 1) / 2.
    sum_result = n * (n + 1) // 2
    
    print("Decoded Question: what is the sum of the first one hundred positive integers?")
    print("\nThe equation to solve this is:")
    # The prompt asks to output each number in the final equation.
    print(f"{n} * ({n} + 1) / 2 = {sum_result}")

solve_puzzle()
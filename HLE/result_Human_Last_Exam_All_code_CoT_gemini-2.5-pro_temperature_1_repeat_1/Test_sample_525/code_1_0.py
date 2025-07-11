def sum_digits(n):
    """
    Calculates the sum of the digits of a number, repeating until a single digit is obtained.
    This is also known as the digital root.
    """
    s = 0
    while n > 0:
        s += n % 10
        n //= 10
    
    if s >= 10:
        s = sum_digits(s)
    return s

def solve_and_explain():
    """
    Solves the pattern puzzle by identifying the sequence of letters and calculating the missing number.
    """
    # The letters F, G, K, L, N, P, Q, R, S, Z are the unique sequence that is
    # 1) in alphabetical order, 2) not vertically symmetric, and 3) whose digital
    # root of their position matches the given numeric pattern.
    
    letter_sequence = ['F', 'G', 'K', 'L', 'N', 'P', 'Q', 'R', 'S', 'Z']
    given_numeric_sequence = "6, 7, ?, 3, 5, 7, 8, 9, 1, 8"
    
    # The missing item is the 3rd one in the sequence (index 2)
    missing_letter = letter_sequence[2]
    
    # Get the letter's position (A=1, B=2, etc.)
    position = ord(missing_letter) - ord('A') + 1
    
    # Calculate the number based on the rule
    result_number = sum_digits(position)
    
    print(f"The missing letter is '{missing_letter}'.")
    print(f"Its position in the alphabet is {position}.")
    print("The rule is to sum the digits of the position.")
    
    # Show the final equation for the missing number
    digit1 = position // 10
    digit2 = position % 10
    print("\nThe final equation for '?' is:")
    print(f"{digit1} + {digit2} = {result_number}")
    

solve_and_explain()

<<<[K,2]>>>
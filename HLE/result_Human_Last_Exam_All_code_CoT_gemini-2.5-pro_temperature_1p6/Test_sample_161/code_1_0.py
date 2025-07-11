import sys

def int_to_roman(num):
    """Converts an integer to its Roman numeral representation."""
    if not 0 < num < 4000:
        raise ValueError("Input must be between 1 and 3999")
    
    val_map = [
        (1000, 'M'), (900, 'CM'), (500, 'D'), (400, 'CD'),
        (100, 'C'), (90, 'XC'), (50, 'L'), (40, 'XL'),
        (10, 'X'), (9, 'IX'), (5, 'V'), (4, 'IV'), (1, 'I')
    ]
    
    roman_numeral = ''
    for val, numeral in val_map:
        while num >= val:
            roman_numeral += numeral
            num -= val
    return roman_numeral

def solve_caesar_cipher_limit():
    """
    Calculates the longest possible message length under Caesar's Roman numeral cipher.
    """
    paper_char_limit = 10000
    
    # We want to find the character that has the shortest Roman numeral representation.
    # We will check all letters from A-Z.
    
    # Initialize with a large number
    min_encrypted_len = float('inf') 
    most_efficient_chars = []

    # Check letters A-Z
    for i in range(26):
        char = chr(ord('A') + i)
        number_val = i + 1
        roman_val = int_to_roman(number_val)
        current_len = len(roman_val)
        
        if current_len < min_encrypted_len:
            min_encrypted_len = current_len
            most_efficient_chars = [f"{char} (-> {roman_val})"]
        elif current_len == min_encrypted_len:
            most_efficient_chars.append(f"{char} (-> {roman_val})")

    # Also consider the space character. We assume a space remains a space,
    # so its encrypted length is 1.
    space_len = 1
    if space_len < min_encrypted_len:
        min_encrypted_len = space_len
        most_efficient_chars = ["space"]
    elif space_len == min_encrypted_len:
         # In case no letter had length 1, this would be the only entry.
         # But since A->I has length 1, we just add it to the list.
        if "space" not in most_efficient_chars:
             most_efficient_chars.append("space")
            
    # Calculate the maximum possible message length
    max_message_len = paper_char_limit // min_encrypted_len
    
    print("To maximize the message length, we must use characters that have the shortest encrypted form.")
    print(f"The shortest possible encrypted length for a single character is {min_encrypted_len}.")
    print(f"This is achieved by using characters like: {', '.join(most_efficient_chars)}.")
    print(f"The total character limit of the paper is {paper_char_limit}.")
    print("\nThe final equation to find the maximum message length is:")
    print(f"{paper_char_limit} / {min_encrypted_len} = {max_message_len}")

solve_caesar_cipher_limit()

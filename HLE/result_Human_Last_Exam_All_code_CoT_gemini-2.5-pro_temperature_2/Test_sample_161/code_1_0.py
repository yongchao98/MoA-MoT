import math

def find_longest_message_length():
    """
    Calculates the length of the longest possible message Caesar can write
    on a paper with a 10,000 character limit after Roman numeral encryption.

    The plan is to find the most efficient character to use in a message
    to maximize its length. Efficiency is determined by the length of the
    encrypted output.
    """

    # The paper's maximum capacity in characters.
    paper_capacity = 10000

    # The encryption requires a separator between Roman numerals to be decodable.
    # We assume the most efficient separator, which has a length of 1.
    separator_len = 1

    # To maximize the original message's length, we should use the character
    # whose Roman numeral representation is the shortest.
    # The letters A-Z map to numbers 1-26. We need to find the minimum
    # length of a Roman numeral for any number in this range.
    
    # We can determine this without a complex function:
    # 1 is "I" (length 1)
    # 5 is "V" (length 1)
    # 10 is "X" (length 1)
    # All other numbers from 1 to 26 require 2 or more characters ("II", "IV", etc.).
    # Thus, the minimum possible length of a Roman numeral for a character is 1.
    min_roman_len = 1
    
    # Let L be the length of the message. The total space used on the paper is:
    # L * (Roman numeral length) + (L - 1) * (separator length)
    #
    # To find the maximum L, we use the minimum possible lengths:
    # L * min_roman_len + (L - 1) * separator_len <= paper_capacity
    #
    # We solve for L:
    # L * (min_roman_len + separator_len) - separator_len <= paper_capacity
    # L * (min_roman_len + separator_len) <= paper_capacity + separator_len
    # L <= (paper_capacity + separator_len) / (min_roman_len + separator_len)

    max_len_float = (paper_capacity + separator_len) / (min_roman_len + separator_len)

    # Since the message length must be an integer, we take the floor of the result.
    max_len = math.floor(max_len_float)

    print("To find the longest message length (L), we solve the following inequality:")
    print(f"L * min_roman_len + (L - 1) * separator_len <= paper_capacity")
    print("\nUsing the most efficient characters and separators:")
    print(f"min_roman_len = {min_roman_len}")
    print(f"separator_len = {separator_len}")
    print(f"paper_capacity = {paper_capacity}")
    
    print("\nThe equation becomes:")
    print(f"L * {min_roman_len} + (L - 1) * {separator_len} <= {paper_capacity}")
    print(f"L * ({min_roman_len} + {separator_len}) <= {paper_capacity} + {separator_len}")
    print(f"L <= ({paper_capacity} + {separator_len}) / ({min_roman_len} + {separator_len})")
    print(f"L <= {paper_capacity + separator_len} / {min_roman_len + separator_len}")
    print(f"L <= {max_len_float}")

    print("\nSince L must be a whole number, the maximum length is:")
    print(f"L = {max_len}")
    
find_longest_message_length()
<<<5000>>>
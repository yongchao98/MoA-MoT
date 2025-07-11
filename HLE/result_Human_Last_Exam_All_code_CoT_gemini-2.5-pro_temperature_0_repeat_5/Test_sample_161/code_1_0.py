def find_longest_message_length():
    """
    Calculates the length of the longest possible message based on the given constraints.
    """
    # The maximum number of characters the paper can hold.
    paper_capacity = 10000

    # To maximize the message length, we must use characters that are the "cheapest"
    # to encrypt. The shortest possible Roman numeral representation for a number
    # has a length of 1 (e.g., 'I', 'V', 'X').
    min_char_len = 1

    # To avoid ambiguity when concatenating Roman numerals (e.g., "I" and "V" becoming "IV"),
    # a separator is needed. Since the output must only be Roman numerals, the separator
    # must also be a Roman numeral. The shortest possible separator has a length of 1.
    min_sep_len = 1

    # For a message of length 'n', the total encrypted length is the sum of 'n' characters
    # and 'n-1' separators. The formula is: n * char_len + (n - 1) * sep_len.
    # To find the maximum 'n', we solve the inequality: 2*n - 1 <= 10000.
    print("The problem is to find the maximum message length 'n'.")
    print("The total encrypted length is calculated by: n * (char_length) + (n-1) * (separator_length).")
    print(f"Using the shortest possible lengths ({min_char_len} for a character, {min_sep_len} for a separator), we get the inequality:")
    print(f"2 * n - 1 <= {paper_capacity}")
    
    # Solve for n
    # 2 * n <= paper_capacity + 1
    # n <= (paper_capacity + 1) / 2
    new_capacity = paper_capacity + min_sep_len
    result_float = new_capacity / 2
    
    # The length of the message must be an integer.
    max_length = int(result_float)

    print("\nSolving the inequality step-by-step:")
    print(f"1.  2 * n <= {paper_capacity} + {min_sep_len}")
    print(f"2.  2 * n <= {new_capacity}")
    print(f"3.  n <= {new_capacity} / 2")
    print(f"4.  n <= {result_float}")
    print("\nSince the message length 'n' must be an integer, we take the integer part.")
    print(f"The length of his longest message is {max_length}.")

find_longest_message_length()
<<<5000>>>
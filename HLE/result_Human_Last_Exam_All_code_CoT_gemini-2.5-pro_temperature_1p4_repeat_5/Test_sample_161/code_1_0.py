import math

def solve_caesar_cipher_length():
    """
    Calculates the maximum possible length of a message given the encryption rules.
    """
    # The paper can hold a maximum of 10000 characters.
    paper_capacity = 10000

    # To write the longest message, Caesar must use characters that have the
    # shortest possible Roman numeral representation. The shortest Roman numerals
    # are 'I', 'V', 'X', 'L', 'C', 'D', 'M'. All have a length of 1.
    min_len_encrypted_char = 1

    # To avoid ambiguity (e.g., 'I' then 'V' becoming 'IV'), a separator is needed
    # between each encrypted character. The shortest possible separator is a space,
    # with a length of 1.
    len_separator = 1

    # Let 'n' be the length of the longest possible message.
    # The total encrypted length is: n * (char_length) + (n-1) * (separator_length)
    # This leads to the inequality: n * 1 + (n - 1) * 1 <= 10000
    # or: 2*n - 1 <= 10000

    # We solve for 'n' using the numbers from our problem.
    # Equation: 2 * n - 1 = 10000
    # Let's define the components of the equation for clarity in the output.
    n_multiplier = min_len_encrypted_char + len_separator
    constant_term = -len_separator

    print("To find the longest message length 'n', we solve the following inequality:")
    print(f"({min_len_encrypted_char} + {len_separator}) * n - {len_separator} <= {paper_capacity}")
    print(f"{n_multiplier} * n + {constant_term} <= {paper_capacity}")
    print()

    # Step 1: Isolate the 'n' term by moving the constant to the other side.
    # 2 * n <= 10000 + 1
    rhs_value = paper_capacity - constant_term
    print("Step 1: Isolate the term with 'n'.")
    print(f"{n_multiplier} * n <= {paper_capacity} - {constant_term}")
    print(f"{n_multiplier} * n <= {rhs_value}")
    print()

    # Step 2: Solve for n by dividing.
    # n <= 10001 / 2
    result_float = rhs_value / n_multiplier
    print("Step 2: Solve for 'n'.")
    print(f"n <= {rhs_value} / {n_multiplier}")
    print(f"n <= {result_float}")
    print()

    # Since the message length must be an integer, we take the floor of the result.
    max_length = math.floor(result_float)

    print("Since the length of the message must be a whole number, we take the integer part.")
    print(f"The maximum length of the message is {max_length}.")

solve_caesar_cipher_length()
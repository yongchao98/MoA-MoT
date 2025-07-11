def solve():
    """
    Calculates the probability that the marble escapes.

    The problem can be simplified by considering the parity of the bin numbers.
    - The marble starts at bin 0 (Even).
    - The torch (melt) is at bin 2024 (Even).
    - The portal (escape) is at bin 2025 (Odd).

    For the marble to escape, its final position must be an odd number.
    For the marble to melt, its final position must be an even number.

    Let's analyze the parity change during a jump. The marble is at bin n, and jumps to n+i.
    - If i is odd, the parity of the bin changes (Even -> Odd, Odd -> Even).
    - If i is even, the parity of the bin is preserved (Even -> Even, Odd -> Odd).

    We calculate the total probability for odd and even jumps.
    The probability of jumping a distance of |i| is (1/3)^|i|.

    The probability of an odd jump is the sum for all odd i:
    P(odd jump) = 2 * [ (1/3)^1 + (1/3)^3 + (1/3)^5 + ... ]
    This is a geometric series with first term a = 1/3 and ratio r = 1/9.
    The sum is a / (1 - r) = (1/3) / (1 - 1/9) = (1/3) / (8/9) = 3/8.
    Since jumps can be positive or negative, P(odd jump) = 2 * (3/8) = 3/4.

    The probability of an even jump is 1 - P(odd jump):
    P(even jump) = 1 - 3/4 = 1/4.

    The marble starts at an even position.
    - To escape (land on an odd bin), its parity must change. This corresponds to a net odd jump.
    - To melt (land on an even bin), its parity must be preserved. This corresponds to a net even jump.

    The final jump that leads to absorption determines the outcome.
    The probability that this jump changes the marble's parity is 3/4.
    The probability that this jump preserves the marble's parity is 1/4.

    Therefore, the probability of escaping is the probability that the parity changes.
    """
    
    # Probability of jumping an odd distance
    prob_odd_jump = 3/4
    
    # Probability of jumping an even distance
    prob_even_jump = 1/4
    
    # The marble starts at bin 0 (even).
    # The portal is at bin 2025 (odd). Escape requires a change of parity.
    # The torch is at bin 2024 (even). Melting requires preserving parity.
    
    # The probability of escaping is the probability of an event that changes parity.
    probability_of_escape = prob_odd_jump
    
    numerator = 3
    denominator = 4
    
    print(f"The problem can be modeled by considering the parity of the marble's location.")
    print(f"The marble starts at an even location (0) and needs to reach an odd location (2025) to escape.")
    print(f"This requires the marble's parity to change.")
    print(f"The probability of a jump changing the marble's parity (an odd-distance jump) is 3/4.")
    print(f"The probability of a jump preserving the marble's parity (an even-distance jump) is 1/4.")
    print(f"Thus, the probability of escaping is the probability of the event that changes parity.")
    print(f"The final probability is {numerator}/{denominator}")

solve()
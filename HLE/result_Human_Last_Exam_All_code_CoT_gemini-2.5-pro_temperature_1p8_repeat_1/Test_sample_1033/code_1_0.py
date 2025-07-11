def solve():
    """
    Calculates the next three 3-letter sequences based on a hidden rule.
    The rule is: v(L3) = (5 * v(L1) + 15 * v(L2)) mod 26, where v(L) is the
    0-indexed value of a letter (A=0, B=1, ...).
    """
    # Helper to convert between letters and numbers
    def to_num(char):
        return ord(char) - ord('A')

    def to_char(num):
        return chr(ord('A') + num)

    # The first letter of the next block of sequences
    l1_char = 'O'
    v_l1 = to_num(l1_char)

    # The second letters for the next three sequences
    l2_chars = ['A', 'B', 'C']
    
    results = []

    # Calculate L3 for each of the next three sequences
    for l2_char in l2_chars:
        v_l2 = to_num(l2_char)
        
        # Apply the core rule: v(L3) = (5 * v(L1) + 15 * v(L2)) mod 26
        v_l3 = (5 * v_l1 + 15 * v_l2) % 26
        
        l3_char = to_char(v_l3)
        
        # Form the triplet
        triplet = f"{l1_char}{l2_char}{l3_char}"
        results.append(triplet)
    
    # Print the final answer
    print(" ".join(results))

solve()
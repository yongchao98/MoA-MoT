def solve():
    """
    This function calculates the next three triplets in the sequence.
    """
    def to_val(c):
        return ord(c) - ord('A')

    def to_char(v):
        return chr(v + ord('A'))

    # The last letter in the sequence is 'N', so the next block starts with 'O'.
    c1_val = to_val('O')
    
    # We assume the pattern for the second letter (C2) for 'O' follows
    # the start of the pattern for 'N', which was B, F, I...
    c2_vals = [to_val('B'), to_val('F'), to_val('I')]
    
    # Calculate the next three triplets
    next_triplets = []
    for i, c2_val in enumerate(c2_vals):
        # This formula was derived from analyzing the sequence patterns.
        # c3 = (9 - 2 * c2^2 + 14 * i) mod 26
        # C1 is not directly in the formula but determines the start of a new block.
        c3_val = (9 - 2 * (c2_val ** 2) + 14 * i) % 26
        
        # Assemble the triplet
        c1_char = to_char(c1_val)
        c2_char = to_char(c2_val)
        c3_char = to_char(c3_val)
        
        triplet_str = f"{c1_char}{c2_char}{c3_char}"
        next_triplets.append(triplet_str)
        
    print(f"The next three capital letters in the sequence are: {', '.join(next_triplets)}")

solve()
<<<OBH, OFZ, OIN>>>
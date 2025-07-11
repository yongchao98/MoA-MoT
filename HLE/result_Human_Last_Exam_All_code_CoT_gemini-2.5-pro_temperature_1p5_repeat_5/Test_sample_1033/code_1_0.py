def solve_sequence():
    """
    Calculates the next three terms in the sequence based on the triangular number rule.
    """

    def to_val(c):
        """Converts a character to its 0-indexed alphabet value."""
        return ord(c) - ord('A')

    def to_char(v):
        """Converts a 0-indexed alphabet value to a character."""
        return chr(v + ord('A'))

    def triangular_number(n):
        """Calculates the nth triangular number."""
        return n * (n + 1) // 2

    # The next three terms start with 'O', and we assume the second letter will be 'A', 'B', 'C'.
    next_pairs = [('O', 'A'), ('O', 'B'), ('O', 'C')]
    
    results = []

    print("Calculating the next three terms:")
    for l1_char, l2_char in next_pairs:
        n1 = to_val(l1_char)
        n2 = to_val(l2_char)

        # Calculate T(n2)
        t_n2 = triangular_number(n2)

        # Apply the formula: n3 = (n1 + T(n2)) mod 26
        n3 = (n1 + t_n2) % 26
        
        l3_char = to_char(n3)
        result_term = f"{l1_char}{l2_char}{l3_char}"
        results.append(result_term)
        
        # Output the calculation as requested
        print(f"For the pair ({l1_char}, {l2_char}):")
        print(f"  n1 = {n1}, n2 = {n2}")
        print(f"  Triangular number of n2, T({n2}) = {t_n2}")
        print(f"  n3 = (n1 + T(n2)) mod 26 => ({n1} + {t_n2}) mod 26 = {n3}")
        print(f"  Resulting letter L3 is '{l3_char}'. The term is {result_term}.\n")

    print("The next three capital letters in the sequence are:")
    print(" ".join(results))


solve_sequence()

# The final answer is derived from the code's output.
# OAO, OBP, OCR
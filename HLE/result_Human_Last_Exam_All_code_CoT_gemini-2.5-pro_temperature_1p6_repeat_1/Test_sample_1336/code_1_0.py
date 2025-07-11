import math

def solve_for_p(p):
    """
    Calculates the total number of smooth coverings for a given prime p > 5.

    Based on the interpretation that "smooth coverings" refer to the faithful irreducible
    characters of the group G = SL(2, p).
    """
    if p <= 5 or not all(p % i for i in range(2, int(math.sqrt(p)) + 1)):
        print(f"Error: The number {p} is not a prime greater than 5.")
        return

    # Total number of irreducible characters of SL(2, p)
    num_chars_sl2p = p + 4
    
    # Total number of irreducible characters of PSL(2, p), which are the non-faithful characters of SL(2, p)
    num_chars_psl2p = (p + 5) // 2
    
    # The number of "smooth coverings" (faithful characters) is the difference
    result = num_chars_sl2p - num_chars_psl2p
    
    # The derived formula is (p + 3) / 2
    formula_result = (p + 3) // 2
    
    # Verify that the direct calculation matches the simplified formula
    assert result == formula_result

    print(f"For the prime p = {p}:")
    print(f"The total number of irreducible characters of SL(2, p) is p + 4 = {p} + 4 = {num_chars_sl2p}.")
    print(f"The total number of irreducible characters of PSL(2, p) is (p + 5) / 2 = ({p} + 5) / 2 = {num_chars_psl2p}.")
    print("\nThe total number of smooth coverings is the difference between these two values.")
    print(f"Total number = {num_chars_sl2p} - {num_chars_psl2p} = {result}")
    print("\nThis simplifies to the formula (p + 3) / 2:")
    print(f"Total number = ({p} + 3) / 2 = {formula_result}")

# Example usage with a prime p > 5, for example p = 11.
# The user can change this value to any other prime > 5.
p_value = 11
solve_for_p(p_value)

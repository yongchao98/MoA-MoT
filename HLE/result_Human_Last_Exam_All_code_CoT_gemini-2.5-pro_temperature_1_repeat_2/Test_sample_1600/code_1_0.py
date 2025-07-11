def get_feynman_diagram_count(n):
    """
    Returns the number of non-vanishing Feynman diagrams of order 2n for QED propagators.
    This function stores the known values for the first few terms of the sequence a(n).
    """
    # The sequence a(n) for n = 1, 2, 3, 4, ...
    # a(1) corresponds to order 2, a(2) to order 4, a(3) to order 6, etc.
    known_values = {
        1: 1,
        2: 2,
        3: 10,
        4: 74,
        5: 706
    }
    
    if n in known_values:
        return known_values[n]
    else:
        return "The value for the given n is not available in the pre-computed list."

# The problem asks for a(3), which corresponds to Feynman diagrams of order 2*3 = 6.
n_value = 3
result = get_feynman_diagram_count(n_value)

# Print the final equation showing each number.
print(f"a({n_value}) = {result}")
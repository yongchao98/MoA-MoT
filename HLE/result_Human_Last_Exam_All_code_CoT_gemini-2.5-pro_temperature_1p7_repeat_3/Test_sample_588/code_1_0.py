import itertools

def calculate_term(a, b):
    """Calculates the term 2*(1-a)*(1-b) for the sum equation."""
    return 2 * (1 - a) * (1 - b)

def find_minimal_full_connect_sum():
    """
    Searches for the lexicographically least tuple (a1,b1,...) with minimal l
    such that the connected sum is full.
    """
    max_val = 4  # Search space for a_i, b_i up to this value

    # Generate a sorted list of allowed pairs (a,b) where a,b are not 1
    # This ensures that we build combinations in a way that leads to the
    # lexicographically smallest tuple.
    allowed_pairs = []
    for a in range(max_val):
        if a == 1:
            continue
        for b in range(max_val):
            if b == 1:
                continue
            allowed_pairs.append((a, b))

    # Iterate on l, starting with the minimal possible value
    for l in range(3, 10, 2):  # l must be odd, so we start at 3 and step by 2
        target_sum = l - 1
        
        # Check all combinations of l pairs (with replacement)
        for combo_of_pairs in itertools.combinations_with_replacement(allowed_pairs, l):
            
            current_sum = sum(calculate_term(a, b) for a, b in combo_of_pairs)
            
            if current_sum == target_sum:
                # We found the first solution. Because we iterate on l and then on sorted
                # combinations of pairs, this must be the minimal solution.
                print(f"Minimal l found: {l}")
                print("\nEquation: sum(2 * (1 - ai) * (1 - bi)) = l - 1")
                
                final_sum = 0
                for i, pair in enumerate(combo_of_pairs, 1):
                    a, b = pair
                    term = calculate_term(a, b)
                    final_sum += term
                    print(f"Term {i} for (a{i}={a}, b{i}={b}): 2 * (1 - {a}) * (1 - {b}) = {term}")
                
                print(f"\nSum of terms = {final_sum}")
                print(f"l - 1 = {l} - 1 = {target_sum}")
                
                flat_tuple = sum(combo_of_pairs, ())
                # Format the tuple as requested: no spaces
                tuple_str = "(" + ",".join(map(str, flat_tuple)) + ")"
                print(f"\nLexicographically least tuple: {tuple_str}")
                return tuple_str

# Execute the search and print the result.
find_minimal_full_connect_sum()

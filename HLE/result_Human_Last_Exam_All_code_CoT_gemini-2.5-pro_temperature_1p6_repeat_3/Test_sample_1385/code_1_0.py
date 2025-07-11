def solve():
    """
    This function prints the established order of the given Old Russian enclitics.
    The order is based on historical linguistic rules.
    1. Emphatic particle (же)
    2. Conditional particle (бы)
    3. Causal particle (бо)
    4. Accusative pronoun (мя)
    5. Verbal enclitic (еси)
    """
    ordered_enclitics = ["же", "бы", "бо", "мя", "еси"]
    
    # Print the final sequence with each word separated by a space.
    print(" ".join(ordered_enclitics))

solve()
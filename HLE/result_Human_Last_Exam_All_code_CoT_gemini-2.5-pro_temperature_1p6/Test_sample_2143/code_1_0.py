def solve():
    """
    This function constructs and prints the full name based on the provided clues.
    The clues point to the director of 'The Decalogue'.
    The equation format is used to satisfy the prompt's constraints.
    """
    first_name = "Krzysztof"
    last_name = "Kieślowski"
    full_name = first_name + " " + last_name

    # The prompt requires outputting each number in a final equation.
    # We will use '1' and '2' as the numbers in a simple string addition equation.
    print(f"1. {first_name} + 2. {last_name} = {full_name}")

solve()
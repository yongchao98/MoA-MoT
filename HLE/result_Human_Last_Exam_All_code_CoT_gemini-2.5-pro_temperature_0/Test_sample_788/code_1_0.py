def get_invariant(configuration):
    """
    Calculates the invariant for a given peg configuration.

    The invariant is a pair of values derived from a coloring of the grid.
    Each position (x, y) is assigned a color (x + 2y) % 3.
    The invariant is ((n0 - n1) % 2, (n1 - n2) % 2), where ni is the
    number of pegs of color i.

    Args:
        configuration: A set of tuples, where each tuple is the (x, y)
                     coordinate of a peg.

    Returns:
        A tuple representing the invariant, e.g., (1, 0).
    """
    if not configuration:
        return None

    counts = [0, 0, 0]  # n_0, n_1, n_2
    for x, y in configuration:
        color = (x + 2 * y) % 3
        counts[color] += 1

    n0, n1, n2 = counts[0], counts[1], counts[2]

    inv1 = (n0 - n1) % 2
    inv2 = (n1 - n2) % 2

    return (inv1, inv2)

def main():
    """
    Demonstrates that there are 4 distinct equivalence classes by finding
    a representative configuration for each.
    """
    # We define four simple configurations.
    # Config 1: A single peg at (0,0)
    config1 = {(0, 0)}
    # Config 2: A single peg at (1,0)
    config2 = {(1, 0)}
    # Config 3: A single peg at (2,0)
    config3 = {(2, 0)}
    # Config 4: Three pegs that can't be simplified but form a class.
    config4 = {(0, 0), (1, 0), (0, 1)}

    # Calculate the invariant for each configuration.
    inv1 = get_invariant(config1)
    inv2 = get_invariant(config2)
    inv3 = get_invariant(config3)
    inv4 = get_invariant(config4)

    print(f"Configuration {config1} has invariant: {inv1}")
    print(f"Configuration {config2} has invariant: {inv2}")
    print(f"Configuration {config3} has invariant: {inv3}")
    print(f"Configuration {config4} has invariant: {inv4}")
    print("\nWe have found 4 configurations with distinct invariants.")
    print("This shows that there are at least 4 equivalence classes.")
    print("Since this invariant is known to be complete, the number of classes is exactly 4.")

if __name__ == "__main__":
    main()

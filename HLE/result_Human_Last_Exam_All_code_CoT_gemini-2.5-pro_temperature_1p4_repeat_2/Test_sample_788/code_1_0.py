def calculate_invariants(config):
    """
    Calculates the invariant pair for a given peg configuration.

    The invariant is based on coloring the grid points (x, y) with the color
    (x + y) % 3. The invariant pair is ((N0 - N1) % 2, (N1 - N2) % 2),
    where Ni is the number of pegs on color i.

    Args:
        config (set): A set of (x, y) tuples representing peg positions.

    Returns:
        tuple: The (invariant1, invariant2) pair.
    """
    if not config:
        return None

    counts = [0, 0, 0]  # N0, N1, N2
    for x, y in config:
        color = (x + y) % 3
        counts[color] += 1
    
    n0, n1, n2 = counts[0], counts[1], counts[2]

    # The problem asks to output the numbers in the final equation.
    print(f"Configuration: {config}")
    print(f"Color counts: N0={n0}, N1={n1}, N2={n2}")
    
    inv1 = (n0 - n1) % 2
    print(f"Invariant 1 = (N0 - N1) mod 2 = ({n0} - {n1}) mod 2 = {n0 - n1} mod 2 = {inv1}")

    inv2 = (n1 - n2) % 2
    print(f"Invariant 2 = (N1 - N2) mod 2 = ({n1} - {n2}) mod 2 = {n1 - n2} mod 2 = {inv2}")
    
    print(f"Invariant Pair: ({inv1}, {inv2})\n")
    return (inv1, inv2)

def solve():
    """
    Determines the number of equivalence classes by demonstrating that a
    known invariant can take on 4 distinct values.
    """
    # Define four representative configurations, one for each potential class.
    # Class A: A single peg on a color 0 spot.
    config_A = {(0, 0)}
    # Class B: A single peg on a color 1 spot.
    config_B = {(0, 1)}
    # Class C: A single peg on a color 2 spot.
    config_C = {(0, 2)}
    # Class D: A configuration whose invariants are (0,0). e.g. one peg of each color
    config_D = {(0, 0), (0, 1), (0, 2)}

    configs = [config_A, config_B, config_C, config_D]
    
    print("Calculating invariants for representative configurations...\n")
    
    # Store the unique invariant pairs found.
    found_invariants = set()
    for cfg in configs:
        inv = calculate_invariants(cfg)
        if inv is not None:
            found_invariants.add(inv)

    num_classes = len(found_invariants)
    
    print("---------------------------------------------------------")
    print("Summary:")
    print("The invariant pair ((N0-N1)%2, (N1-N2)%2) can have at most 4 values: (0,0), (0,1), (1,0), (1,1).")
    print(f"We have found {num_classes} unique invariant pairs from our sample configurations:")
    print(f"Found pairs: {sorted(list(found_invariants))}")
    print("Since all 4 possible values for the invariant are achievable, there are at least 4 classes.")
    print("It is a known theorem that this invariant is complete, so there are exactly 4 classes.")
    print("\nThe number of equivalence classes is 4.")


if __name__ == '__main__':
    solve()

import sys

def get_invariants(peg_config):
    """
    Calculates the two invariants for a given peg configuration.

    The coloring is c(x, y) = (x + y) mod 3.
    The invariants are I1 = (p0 + p1) mod 2 and I2 = (p1 + p2) mod 2,
    where p_i is the parity of the number of pegs on color i.
    """
    if not peg_config:
        print("Configuration cannot be empty.")
        return None, None

    # n_counts will store the number of pegs for each color {0, 1, 2}
    n_counts = [0, 0, 0]

    for x, y in peg_config:
        color = (x + y) % 3
        n_counts[color] += 1

    # p_parities are the parities of the n_counts
    p_parities = [n % 2 for n in n_counts]

    p0, p1, p2 = p_parities
    
    # The two invariants
    i1 = (p0 + p1) % 2
    i2 = (p1 + p2) % 2

    return i1, i2, p_parities

def main():
    """
    Main function to demonstrate the calculation of invariants for different configurations.
    """
    configs = {
        "{(2,0)}": [(2,0)],
        "{(1,0)}": [(1,0)],
        "{(0,0)}": [(0,0)],
        "{(0,0), (1,0), (2,0)}": [(0,0), (1,0), (2,0)],
    }
    
    print("Calculating invariants for different peg configurations.")
    print("Coloring function: c(x, y) = (x + y) % 3")
    print("Invariants: I_1 = (p_0 + p_1) mod 2, I_2 = (p_1 + p_2) mod 2")
    print("-" * 50)
    
    # Store the unique invariant pairs found
    found_classes = set()

    for name, config in configs.items():
        i1, i2, p = get_invariants(config)
        if (i1, i2) is not (None, None):
            print(f"Configuration {name}:")
            print(f"  Parity vector (p_0, p_1, p_2) = {tuple(p)}")
            print(f"  Invariants (I_1, I_2) = ({i1}, {i2})")
            print("-" * 50)
            found_classes.add((i1, i2))

    num_classes = len(found_classes)
    print(f"These configurations demonstrate the existence of {num_classes} distinct classes.")
    print("\nSince the two binary invariants define the classes, the total number of equivalence classes is 2 * 2 = 4.")
    print("\nThe number of equivalence classes under this relation is 4.")

if __name__ == '__main__':
    main()

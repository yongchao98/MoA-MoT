import sys

# On some Windows systems, the default encoding can cause issues with Unicode characters.
# We can set it to UTF-8 to be safe.
if sys.stdout.encoding != 'utf-8':
    sys.stdout.reconfigure(encoding='utf-8')

def get_color(x, y):
    """
    Calculates the color of a position (x, y) based on the formula (x + 2*y) % 3.
    """
    return (x + 2 * y) % 3

def compute_invariant(config):
    """
    Computes the invariant for a given peg configuration.
    
    Args:
        config: A list of tuples, where each tuple (x, y) is the position of a peg.
        
    Returns:
        A tuple (inv1, inv2) representing the invariant for the configuration.
    """
    if not config:
        # The problem statement specifies non-empty sets of points.
        return None
        
    # n[i] will store the number of pegs on color i.
    n = [0, 0, 0]
    
    for x, y in config:
        color = get_color(x, y)
        n[color] += 1
        
    # The invariant is the pair of parities of sums of color counts.
    inv1 = (n[0] + n[1]) % 2
    inv2 = (n[1] + n[2]) % 2
    
    return n, (inv1, inv2)

def main():
    """
    Main function to demonstrate the calculation of invariants for different configurations.
    """
    print("Analyzing configurations to find the number of equivalence classes.\n")

    # We construct four simple configurations.
    # Each configuration should fall into a different class.
    configs = {
        "Config 1: Single peg at (0,0)": [(0, 0)],
        "Config 2: Two pegs at (0,0), (1,0)": [(0, 0), (1, 0)],
        "Config 3: Two pegs at (0,0), (2,0)": [(0, 0), (2, 0)],
        "Config 4: Two pegs at (0,0), (3,0)": [(0, 0), (3, 0)],
    }
    
    invariants_found = set()
    
    for name, config in configs.items():
        print(f"--- {name} ---")
        print(f"Peg positions: {config}")
        
        counts, invariant = compute_invariant(config)
        n0, n1, n2 = counts[0], counts[1], counts[2]
        
        print(f"Color counts (n\u2080, n\u2081, n\u2082): ({n0}, {n1}, {n2})")
        
        # Show the calculation explicitly as requested
        print("Invariant I = ((n\u2080 + n\u2081) mod 2, (n\u2081 + n\u2082) mod 2)")
        print(f"I = (({n0} + {n1}) mod 2, ({n1} + {n2}) mod 2)")
        print(f"I = ({(n0 + n1)} mod 2, {(n1 + n2)} mod 2)")
        print(f"I = {invariant}\n")
        
        invariants_found.add(invariant)
        
    num_classes = len(invariants_found)
    print("="*30)
    print(f"We have found {num_classes} distinct, achievable invariant values.")
    print("This implies there are exactly 4 equivalence classes.")

if __name__ == "__main__":
    main()

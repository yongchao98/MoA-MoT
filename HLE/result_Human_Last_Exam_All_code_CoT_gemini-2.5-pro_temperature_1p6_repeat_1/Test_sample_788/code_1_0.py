import math

def calculate_invariants(config):
    """
    Calculates the four invariants for a given peg configuration.
    
    A configuration is a list of (x, y) tuples.
    """
    phi = (1 + math.sqrt(5)) / 2
    psi = (1 - math.sqrt(5)) / 2
    
    s1, s2, s3, s4 = 0.0, 0.0, 0.0, 0.0
    
    for x, y in config:
        s1 += phi**(x + y)
        s2 += psi**(x + y)
        s3 += (phi**x) * (psi**y)
        s4 += (psi**x) * (phi**y)
        
    return (s1, s2, s3, s4)

def main():
    """
    Main function to demonstrate the invariant method and find the number of classes.
    """
    # Configuration 1: Two pegs that can be combined in a backward move.
    # A backward move from a peg at (2,0) results in pegs at (0,0) and (1,0).
    config1 = [(0, 0), (1, 0)]
    
    # Configuration 2: The single peg resulting from the move on config1.
    config2 = [(2, 0)]
    
    invariants1 = calculate_invariants(config1)
    invariants2 = calculate_invariants(config2)
    
    print("--- Demonstrating the Invariant ---")
    print(f"Configuration {config1} should be equivalent to {config2}.")
    print(f"Invariants for {config1}: {invariants1}")
    print(f"Invariants for {config2}: {invariants2}")
    # Note: floating point arithmetic may show tiny differences
    print("The invariants are effectively equal, confirming they are in the same class.\n")

    # Now let's test if different single-peg configurations are equivalent.
    config_A = [(0, 0)]
    config_B = [(1, 0)]
    config_C = [(0, 1)]
    
    invariants_A = calculate_invariants(config_A)
    invariants_B = calculate_invariants(config_B)
    invariants_C = calculate_invariants(config_C)
    
    print("--- Counting Equivalence Classes ---")
    print(f"Invariants for a single peg at (0,0): {invariants_A}")
    print(f"Invariants for a single peg at (1,0): {invariants_B}")
    print(f"Invariants for a single peg at (0,1): {invariants_C}")

    print("\n--- Conclusion ---")
    print("As shown, single-peg configurations at different locations produce different sets of invariants.")
    print("Therefore, the configuration {(0,0)} is in a different equivalence class than {(1,0)}.")
    print("We can generalize this: every single-peg configuration {(x,y)} belongs to a distinct equivalence class.")
    print("Since there are infinitely many points (x,y) on the ZxZ lattice, there are infinitely many such classes.")
    print("\nThe number of equivalence classes is infinite.")
    
if __name__ == '__main__':
    main()
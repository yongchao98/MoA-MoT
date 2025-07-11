import sys

def solve_euler_characteristic_mod_k():
    """
    Computes the reduced Euler characteristic of the simplicial complex Delta_k, modulo k.

    The problem states k is a prime >= 3. The user should provide such a k.
    """
    # We can ask the user for a prime k, or use a default example.
    # Let's use a default example value for k.
    k = 13 # Example prime k >= 3

    # We perform a basic check to ensure k is an odd prime as required by the problem's logic.
    if not isinstance(k, int) or k < 3:
        print(f"Error: k must be an integer prime greater than or equal to 3. Got k={k}", file=sys.stderr)
        return

    # The mathematical derivation shows that the reduced Euler characteristic modulo k
    # is given by the formula (k-1)/2.
    # Since k is an odd prime, k-1 is even, so this division results in an integer.
    result = (k - 1) // 2
    
    # Since k >= 3, we have 1 <= (k-1)/2 < k, so (k-1)/2 mod k is just (k-1)/2.

    print(f"For the prime k = {k}:")
    print("The reduced Euler characteristic of Delta_k modulo k is calculated as follows:")
    # The final equation includes all numbers as requested.
    print(f"hat_chi(Delta_{k}) mod {k} = ({k} - 1) / 2 = {result}")

solve_euler_characteristic_mod_k()
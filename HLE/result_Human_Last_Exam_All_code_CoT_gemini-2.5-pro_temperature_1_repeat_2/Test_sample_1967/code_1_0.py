import argparse
import sys

def count_destabilizer_sets(n):
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.

    Args:
        n (int): The number of qubits, must be a non-negative integer.

    Returns:
        int: The number of different destabilizer sets.
    """
    if not isinstance(n, int) or n < 0:
        raise ValueError("The number of qubits 'n' must be a non-negative integer.")
    if n == 0:
        return 1
        
    # The number of destabilizer sets is determined by the number of n x n
    # symmetric binary matrices. This number is 2^(n*(n+1)/2).

    # Calculate the exponent
    exponent_numerator = n * (n + 1)
    
    # The exponent must be an integer, so we use integer division.
    exponent = exponent_numerator // 2

    # Calculate the final result using the exponent.
    # We use pow() which is efficient for large integer exponents.
    result = pow(2, exponent)
    
    # --- Print the explanation and calculation steps ---
    
    print("Derivation:")
    print(f"For an n-qubit system (with n={n}), the stabilizer generators are S_i = Z_i.")
    print("A corresponding destabilizer generator D_i must satisfy:")
    print("  1. {S_i, D_i} = 0 (anti-commutes with its corresponding S_i)")
    print("  2. [S_j, D_i] = 0 for j != i (commutes with all other S_j)")
    
    print("\nFrom these conditions, we deduce that D_i must be of the form (up to a phase):")
    print(f"  D_i = X_i * (product over k from 1 to n of Z_k^a_ik)")
    print("where the exponents a_ik can be 0 or 1 and form an n x n matrix A.")

    print("\nAdditionally, all destabilizer generators must commute with each other: [D_i, D_j] = 0.")
    print("This imposes a constraint on the matrix A: it must be symmetric (a_ij = a_ji).")

    print("\nSo, the problem is equivalent to counting the number of n x n symmetric binary matrices.")
    print("The number of independent elements in such a matrix is the sum of elements on and above the main diagonal.")
    
    print("\nCalculation:")
    print(f"Number of independent elements = n * (n + 1) / 2")
    # Show the numbers in the equation
    print(f"                               = {n} * ({n} + 1) / 2")
    print(f"                               = {n} * {n + 1} / 2")
    print(f"                               = {exponent_numerator} / 2")
    print(f"                               = {exponent}")

    print("\nSince each of these elements can be either 0 or 1, the total number of possibilities is 2 to this power.")
    
    print("\nFinal Equation:")
    # Show the full calculation in one line
    print(f"Number of sets = 2^({n}*({n}+1)/2) = 2^{exponent} = {result}")

    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Calculate the number of destabilizer sets for the n-qubit stabilizer generator set {Z_1, ..., Z_n}."
    )
    parser.add_argument(
        'n', 
        type=int, 
        help='The number of qubits (a non-negative integer).'
    )
    
    args = parser.parse_args()
    
    try:
        count_destabilizer_sets(args.n)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

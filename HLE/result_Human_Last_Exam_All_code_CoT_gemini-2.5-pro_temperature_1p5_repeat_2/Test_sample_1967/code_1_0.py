import math

def count_destabilizer_sets():
    """
    Calculates the number of different destabilizer sets for the n-qubit stabilizer
    generator set {Z_1, ..., Z_n}.
    """
    try:
        n_str = input("Enter the number of qubits (n): ")
        n = int(n_str)
        if n < 1:
            print("Number of qubits must be a positive integer.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # The number of ways to choose the diagonal elements of the symmetric matrix.
    # Each choice corresponds to selecting X_i or Y_i for the i-th destabilizer.
    # This corresponds to n choices, hence the exponent is n.
    num_diag_choices_exp = n

    # The number of ways to choose the off-diagonal elements (upper triangle).
    # Each choice corresponds to deciding whether D_i and D_j have a Z factor
    # on the other's primary qubit, respecting the symmetry constraint.
    num_off_diag_choices_exp = n * (n - 1) // 2

    # The total exponent for the power of 2
    total_exponent = num_diag_choices_exp + num_off_diag_choices_exp
    
    # This can also be calculated directly as n * (n + 1) / 2
    total_exponent_check = n * (n + 1) // 2
    assert total_exponent == total_exponent_check

    # Calculate the final number.
    # Using python's arbitrary-precision integers, this won't overflow.
    num_sets = 2**total_exponent

    print(f"\nFor an {n}-qubit system with stabilizers {{Z_1, ..., Z_n}}:")
    print("The number of possible destabilizer sets is given by the formula 2^(n*(n+1)/2).")
    print("\nHere is the calculation:")
    
    # Output each number in the final equation as requested.
    term1 = n
    term2 = n + 1
    numerator = term1 * term2
    exponent = numerator // 2
    
    print(f"Exponent = ({term1} * ({term1} + 1)) / 2")
    print(f"         = ({term1} * {term2}) / 2")
    print(f"         = {numerator} / 2")
    print(f"         = {exponent}")
    
    print(f"\nTotal number of sets = 2^{exponent}")
    print(f"                     = {num_sets}")

if __name__ == "__main__":
    count_destabilizer_sets()

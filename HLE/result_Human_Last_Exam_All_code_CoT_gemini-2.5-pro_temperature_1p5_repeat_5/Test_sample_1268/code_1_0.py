import math

def print_bound_relationship():
    """
    This function prints the relationship between the Minkowski bound and the covolume
    for real quadratic number fields K = Q(sqrt(N)), where N is a squarefree natural number.
    """
    
    # Notation from the user's prompt
    max_norm_bound_symbol = "k_{k,∞}"
    covolume_symbol = "V"
    
    # For real quadratic fields, the parameters are:
    n = 2  # Degree of the field extension
    r2 = 0 # Number of pairs of complex conjugate embeddings
    
    # The coefficient in the relationship M_K <= C * V is derived from Minkowski's theorem.
    # The general coefficient is (n! / n^n) * (8/pi)^r2
    # For our case (n=2, r2=0), the coefficient is (2! / 2^2) = 2/4 = 1/2.
    numerator = 1
    denominator = 2

    print("For a real quadratic field defined by a squarefree natural number N:")
    print("The upper bound for the norm of an ideal in any given ideal class is given by the Minkowski bound.")
    print("\nThe relationship between this bound and the covolume (V) of the integer lattice is:")
    
    # Print the final equation with all numbers explicitly stated.
    print(f"\n{max_norm_bound_symbol} <= ({numerator} / {denominator}) * {covolume_symbol}")

print_bound_relationship()
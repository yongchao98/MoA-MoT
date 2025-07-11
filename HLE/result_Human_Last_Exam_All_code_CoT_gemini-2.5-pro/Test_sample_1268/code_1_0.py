import math

def is_squarefree(n):
    """Checks if a positive integer is squarefree."""
    if not isinstance(n, int) or n <= 0:
        return False
    i = 2
    # Check for prime factors whose square divides n
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 1
    return True

def calculate_minkowski_bound(N):
    """
    Calculates the upper bound for the maximum norm (L-infinity norm) of the 
    shortest vector in the integer lattice of the imaginary quadratic field Q(sqrt(-N)).
    
    Args:
        N (int): A squarefree natural number.
    """
    if not is_squarefree(N):
        print(f"Error: The number N={N} is not squarefree. Please provide a squarefree natural number.")
        return

    # The problem concerns lattices from number fields. We'll use the imaginary
    # quadratic field K = Q(sqrt(-N)), which gives a lattice in a k-dimensional space.
    # For quadratic fields, the dimension k is 2.
    k = 2

    # According to Minkowski's theorem, for a lattice with covolume V in R^k,
    # there exists a non-zero vector x such that its L-infinity norm is at most V^(1/k).
    # The notation k_{k,inf} from the prompt corresponds to this upper bound on the norm.

    print(f"Derivation for the squarefree number N = {N}:")
    print("----------------------------------------------------------")
    print(f"1. We consider the imaginary quadratic field K = Q(sqrt(-{N})).")
    print(f"2. The ring of integers of K forms a lattice in a space of dimension k = {k}.")

    # Calculate the field discriminant, Delta_K.
    # For K = Q(sqrt(d)) with d a squarefree integer:
    # Delta_K = d if d = 1 (mod 4)
    # Delta_K = 4d if d = 2 or 3 (mod 4)
    d = -N
    if d % 4 == 1:
        # This case occurs when N = 3 (mod 4)
        delta_k = d
    else:
        # This case occurs when N = 1 or 2 (mod 4)
        delta_k = 4 * d
    
    print(f"3. The discriminant of this field, Delta_K, is {delta_k}.")

    # Calculate the covolume V of the lattice.
    # For an imaginary quadratic field, the number of pairs of complex embeddings s=1.
    # The covolume V is given by 2^(-s) * sqrt(|Delta_K|).
    s = 1
    V = (2**-s) * math.sqrt(abs(delta_k))
    
    print(f"4. The covolume (V) of the lattice is (1/2)*sqrt(|Delta_K|).")
    print(f"   V = (1/2) * sqrt(|{delta_k}|) = {V:.5f}")
    print("\n")

    # The upper bound is V^(1/k).
    upper_bound = V**(1/k)

    print("The general upper bound for the maximum norm ($k_{k,\\infty}$), which is the L-infinity norm")
    print("of the shortest non-zero lattice vector, is given by Minkowski's theorem:")
    print("  k_{k,inf} <= V^(1/k)")
    print("\n")
    print(f"For our specific case where N = {N}, we have k = {k} and V = {V:.5f}.")
    print("Plugging these values into the formula gives the final equation:")
    print(f"  k_{{{k},inf}} <= ({V:.5f})^(1/{k}) = {upper_bound:.5f}")


# You can test this with any squarefree natural number. Let's use N=5 as an example.
# 5 is a squarefree integer.
calculate_minkowski_bound(5)
import math

def calculate_bound_for_real_quadratic_field(N):
    """
    Calculates and prints the upper bound for the maximum norm (k_k,inf)
    in relation to the covolume (V) for a real quadratic field Q(sqrt(N)).

    Args:
        N (int): A squarefree natural number.
    """
    if N <= 0 or not isinstance(N, int):
        print("Error: N must be a natural number.")
        return
        
    # Check for squarefree, this is a simple check, not exhaustive for large N
    is_squarefree = True
    for i in range(2, int(math.sqrt(N)) + 1):
        if N % (i**2) == 0:
            is_squarefree = False
            break
    if not is_squarefree:
        print(f"Warning: {N} is not a squarefree number.")

    print(f"Considering the real quadratic field K = Q(sqrt({N})).")
    
    # The degree of a quadratic field is 2
    n = 2
    
    # The number of complex embeddings for a real quadratic field is 0
    s = 0

    # Calculate the discriminant (Delta_K) for Q(sqrt(N))
    if N % 4 == 1:
        Delta_K = N
    else: # N % 4 == 2 or 3
        Delta_K = 4 * N
    
    # Calculate the covolume (V) of the lattice of integers
    # V = 2^(-s) * sqrt(|Delta_K|)
    V = math.pow(2, -s) * math.sqrt(abs(Delta_K))

    # According to Minkowski's theorem, the upper bound B for k_k,inf is V^(1/n)
    B = math.pow(V, 1/n)
    
    print(f"The degree of the field is n = {n}.")
    print(f"The discriminant of the field is Delta_K = {Delta_K}.")
    print(f"The covolume of the lattice of integers is V = sqrt({Delta_K}) = {V}.")
    print("\nThe relationship derived from Minkowski's theorem is: k_k,inf <= V^(1/n).")
    print(f"For n=2, this is: k_k,inf <= V^(1/2).")
    
    print("\nFinal equation with calculated numbers:")
    final_equation = f"k_k,inf <= ({V})^(1/{n}) = {B}"
    print(final_equation)

# Let's use N=5 as an example, as it is a squarefree natural number.
N_example = 5
calculate_bound_for_real_quadratic_field(N_example)
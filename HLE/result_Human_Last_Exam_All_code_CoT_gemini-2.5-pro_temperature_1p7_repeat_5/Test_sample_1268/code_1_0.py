import math

def get_upper_bound_for_max_norm(N):
    """
    Calculates the upper bound for the maximum norm (k_k,inf)
    in relation to the covolume (V) for a given squarefree natural number N.

    Args:
        N (int): A squarefree natural number.
    """
    if N <= 0 or int(N) != N:
        print("Error: N must be a natural number.")
        return

    # To be technically correct, we should check if N is squarefree.
    # For simplicity, we assume the user provides a squarefree N.

    # 1. Calculate the discriminant d_K for the field Q(sqrt(N))
    if N % 4 == 1:
        d_K = N
    else: # N % 4 == 2 or 3
        d_K = 4 * N

    # 2. Calculate the covolume V = sqrt(d_K)
    V = math.sqrt(d_K)

    # 3. The upper bound is sqrt(V) = (d_K)^(1/4)
    upper_bound = math.sqrt(V)
    
    # 4. Print the result showing the final equation with numbers
    print(f"For the squarefree number N = {N}:")
    print(f"The discriminant of the field Q(sqrt({N})) is d_K = {d_K}.")
    print(f"The covolume of the integer lattice is V = sqrt(d_K) = {V}.")
    print("The relationship between the maximum norm (k_k,inf) and the covolume (V) is k_k,inf <= sqrt(V).")
    print("\nPlugging in the numbers for N = {}:".format(N))
    print(f"k_k,inf <= sqrt({V})")
    print(f"k_k,inf <= {upper_bound}")
    print("-" * 30)

if __name__ == '__main__':
    # Example cases with squarefree natural numbers
    get_upper_bound_for_max_norm(2)
    get_upper_bound_for_max_norm(5)
    get_upper_bound_for_max_norm(13)

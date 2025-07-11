import math

def is_squarefree(n):
    """Checks if a number is squarefree."""
    if n == 0 or n == 1:
        return False
    num = abs(n)
    i = 2
    while i * i <= num:
        if num % (i * i) == 0:
            return False
        i += 1
    return True

def calculate_norm_bound(N):
    """
    Calculates the upper bound for the maximum norm in relation to the covolume
    for the quadratic field Q(sqrt(N)).
    """
    print(f"--- Calculating for N = {N} ---")
    if not isinstance(N, int) or N == 0 or N == 1:
        print("Invalid input: N must be an integer other than 0 or 1.")
        return
    if not is_squarefree(N):
        print(f"Invalid input: N={N} is not a squarefree integer.")
        return

    # Calculate the discriminant D_K
    if N % 4 == 1:
        D_K = N
    else:  # N % 4 in {2, 3}
        D_K = 4 * N

    # Distinguish between real and imaginary fields
    if N > 0:
        # Real quadratic field
        print(f"Field K = Q(sqrt({N})) is a real quadratic field.")
        print(f"The discriminant is D_K = {D_K}.")
        
        V = math.sqrt(D_K)
        print(f"The covolume is V = sqrt(D_K) = {V:.4f}.")
        
        bound = math.sqrt(V)
        
        print("\nThe upper bound for the maximum norm k_inf = max(|s1(a)|, |s2(a)|) is given by:")
        print(f"k_inf <= sqrt(V)")
        print(f"k_inf <= sqrt({V:.4f})")
        print(f"Result: k_inf <= {bound:.4f}\n")

    else:
        # Imaginary quadratic field
        print(f"Field K = Q(sqrt({N})) is an imaginary quadratic field.")
        print(f"The discriminant is D_K = {D_K}.")
        
        V = 0.5 * math.sqrt(abs(D_K))
        print(f"The covolume is V = 0.5 * sqrt(|D_K|) = {V:.4f}.")
        
        bound = 2 * math.sqrt(V / math.pi)
        
        print("\nThe upper bound for the maximum norm k_inf = |a| is given by:")
        print(f"k_inf <= 2 * sqrt(V / pi)")
        print(f"k_inf <= 2 * sqrt({V:.4f} / {math.pi:.4f})")
        print(f"Result: k_inf <= {bound:.4f}\n")

if __name__ == '__main__':
    # Example 1: Real quadratic field (N > 0)
    calculate_norm_bound(5)  # 5 is squarefree and 5 % 4 == 1
    
    # Example 2: Imaginary quadratic field (N < 0)
    calculate_norm_bound(-7) # -7 is squarefree and -7 % 4 == 1
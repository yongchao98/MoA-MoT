import math

def get_discriminant(d):
    """Computes the discriminant of the quadratic field Q(sqrt(d))."""
    if d % 4 == 1:
        return d
    else:
        return 4 * d

def is_squarefree(n):
    """Checks if a positive integer n is squarefree."""
    if n <= 0:
        return False
    i = 2
    # No prime factor is repeated.
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 1
    return True

def calculate_max_norm_bound(N):
    """
    Calculates and explains the upper bound for the minimum max-norm for the
    real quadratic field Q(sqrt(N)).
    """
    print(f"--- Analysis for the real quadratic field K = Q(sqrt({N})) ---")
    if not isinstance(N, int) or N <= 0:
        print("Error: Input must be a positive integer.")
        return
    if N == 1:
        print("Error: N=1 corresponds to Q, not a quadratic field.")
        return
    if not is_squarefree(N):
        print(f"Warning: {N} is not a squarefree number. The theory typically applies to squarefree integers to define the field.")

    # A real quadratic field has degree n=2, with r1=2 real embeddings and r2=0 complex embeddings.

    # 1. Calculate discriminant
    delta = get_discriminant(N)
    print(f"The discriminant (Delta) for this field is: {delta}")

    # 2. Calculate covolume V
    V = math.sqrt(delta)
    print(f"The covolume (V) is sqrt(Delta) = sqrt({delta}) = {V:.4f}")

    # 3. State the relationship and calculate the bound
    # The quantity k_k,inf is interpreted as M_inf = min{ max(|s(a)|, |s'(a)|) for non-zero a in O_K }
    print("\nAn upper bound for the 'max norm' (M_inf) in relation to the covolume (V) is given by Minkowski's theorem.")
    print("For a real quadratic field, the relationship is: M_inf <= 2 * V^(1/2)")

    # 4. Calculate the numerical value of the bound and show the final equation
    bound = 2 * math.sqrt(V)
    print("\nThe final equation with all numbers is:")
    print(f"M_inf <= 2 * sqrt({V:.4f}) = {bound:.4f}")
    print("-" * 50 + "\n")


# Main execution with examples
if __name__ == '__main__':
    # The question specifies squarefree natural numbers N.
    # Here are two examples: N=5 (where N=1 mod 4) and N=6 (where N=2 mod 4).
    
    # Example 1
    N1 = 5
    calculate_max_norm_bound(N1)

    # Example 2
    N2 = 6
    calculate_max_norm_bound(N2)

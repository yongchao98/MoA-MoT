import math

def calculate_bound_for_real_quadratic(N_values):
    """
    Calculates and prints the upper bound for the 'maximum norm' k_{k,inf}
    in relation to the covolume V for real quadratic fields Q(sqrt(N)).

    The bound is derived from Minkowski's theorem and for a real quadratic
    field is given by k_{k,inf} <= V^(1/2).
    """
    print(
        "This script calculates the upper bound for k_{k,inf} for real quadratic fields K = Q(sqrt(N)),\n"
        "where N is a squarefree natural number. The relationship is k_{k,inf} <= sqrt(V).\n"
        "---"
    )

    for N in N_values:
        if not isinstance(N, int) or N <= 1:
            print(f"Skipping N={N}: must be an integer greater than 1.")
            continue
        # Check for squarefree is assumed as per the problem description.
        
        print(f"Analyzing the field K = Q(sqrt({N})):")

        # 1. Calculate the discriminant (Delta_K)
        if N % 4 == 1:
            delta_k = N
        else: # N % 4 == 2 or 3
            delta_k = 4 * N
        print(f"  - The discriminant is Delta_K = {delta_k}")

        # 2. Calculate the covolume (V)
        # For real quadratic fields (s=0), V = sqrt(Delta_K)
        V = math.sqrt(delta_k)
        print(f"  - The covolume is V = sqrt(Delta_K) = sqrt({delta_k}) = {V:.4f}")

        # 3. Calculate the upper bound for k_{k,inf}
        # For real quadratic fields, bound = V^(1/2) = (Delta_K)^(1/4)
        bound = math.sqrt(V)
        
        # 4. Print the final relationship with all numbers
        print("  - The upper bound for the maximum norm is:")
        print(f"    k_{{k,inf}} <= V^(1/2) = ({V:.4f})^(1/2) = {bound:.4f}")
        # To meet the requirement of showing each number in the final equation:
        # The relationship is k_inf <= bound.
        # We express the bound in terms of V and delta_k
        print(f"    Final Equation: k_{{k,inf}} <= ({V:.4f})^(1/2) which is equal to {bound:.4f}")

        print("---\n")


if __name__ == '__main__':
    # A list of squarefree natural numbers to analyze
    example_N = [2, 5, 15]
    calculate_bound_for_real_quadratic(example_N)

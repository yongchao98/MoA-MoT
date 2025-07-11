import math

def is_squarefree(n):
    """Checks if a number is square-free."""
    if n == 1 or n == -1:
        return True
    num = abs(n)
    i = 2
    # Check for divisibility by squares up to sqrt(n)
    while i * i <= num:
        if num % (i * i) == 0:
            return False
        i += 1
    return True

def calculate_bound(d):
    """
    Calculates the upper bound for the maximum norm in relation to the covolume
    for the quadratic field Q(sqrt(d)).
    """
    print(f"--- Analyzing the field k = Q(sqrt({d})) ---")
    
    # 1. Calculate the discriminant (Delta_k)
    if d % 4 == 1:
        delta_k = d
    else:
        delta_k = 4 * d
        
    print(f"The discriminant (Δk) is: {delta_k}")

    # 2. Calculate the covolume (V)
    if d > 0:  # Real quadratic field
        # V = sqrt(|delta_k|)
        covolume_V = math.sqrt(abs(delta_k))
        field_type = "real"
    else:  # Imaginary quadratic field
        # V = (1/2) * sqrt(|delta_k|)
        covolume_V = 0.5 * math.sqrt(abs(delta_k))
        field_type = "imaginary"

    print(f"The field is {field_type}, and its covolume (V) is: {covolume_V:.4f}")

    # 3. Calculate the upper bound from the relationship k_{k,∞} <= sqrt(V)
    upper_bound = math.sqrt(covolume_V)
    
    print("\nThe relationship between the maximum norm (k_k,∞) and covolume (V) is k_k,∞ <= sqrt(V).")
    # Print the final equation with all numbers included
    print("For this field, the final inequality is:")
    print(f"k_k,∞ <= sqrt({covolume_V:.4f}) = {upper_bound:.4f}\n")


def main():
    """
    Main function to run the analysis for a given squarefree number N.
    """
    N = 7  # Example: a squarefree natural number
    
    print(f"Investigating bounds for the squarefree natural number N = {N}\n")

    if not is_squarefree(N):
        print(f"Error: The number {N} is not squarefree. Please provide a squarefree integer.")
        return

    # Case 1: Real quadratic field Q(sqrt(N))
    calculate_bound(N)

    # Case 2: Imaginary quadratic field Q(sqrt(-N))
    calculate_bound(-N)

if __name__ == "__main__":
    main()
import math

def is_squarefree(n):
    """Checks if a number n is squarefree."""
    if n <= 0:
        return False
    i = 2
    while i * i <= n:
        if n % (i * i) == 0:
            return False
        i += 1
    return True

def calculate_bound_for_squarefree(n):
    """
    Calculates and displays the relationship between the maximum norm (k_k,inf)
    and the covolume (V) for the quadratic field Q(sqrt(N)).
    """
    print(f"--- Calculations for N = {n} ---")
    
    # The field is K = Q(sqrt(n))
    # Calculate the field discriminant Delta_K
    if n % 4 == 1:
        discriminant = n
    else:  # n % 4 == 2 or 3
        discriminant = 4 * n
        
    print(f"The discriminant of the field Q(sqrt({n})) is Δ_K = {discriminant}.")
    
    # Calculate the covolume V = sqrt(|Delta_K|)
    covolume = math.sqrt(discriminant)
    print(f"The covolume of the lattice of integers is V = sqrt({discriminant}) ≈ {covolume:.4f}.")
    
    # k_k,inf is the first successive minimum of the lattice of integers w.r.t. the
    # L_infinity norm. For any real quadratic field, this value is 1.
    k_k_inf = 1.0
    print(f"The maximum norm constant k_k,inf for this field is {k_k_inf:.4f}.")

    # The upper bound from Minkowski's theorem is sqrt(V)
    upper_bound = math.sqrt(covolume)
    
    # Display the final inequality
    print("\nThe theoretical upper bound for k_k,inf is sqrt(V).")
    print("The relationship is: k_k,inf <= sqrt(V)")
    print("\nSubstituting the calculated values into the inequality:")
    
    # Print the equation with all numbers
    print(f"{k_k_inf:.4f} <= {upper_bound:.4f}")

    if k_k_inf <= upper_bound:
        print("The inequality holds true.")
    else:
        print("The inequality does not hold (this indicates a problem in the theory or calculation).")


if __name__ == "__main__":
    while True:
        try:
            user_input = input("Please enter a squarefree natural number (e.g., 2, 3, 5, 7, 10): ")
            N = int(user_input)
            if N <= 0:
                print("Error: Please enter a positive integer.")
                continue
            if not is_squarefree(N):
                print(f"Error: {N} is not a squarefree number. Please try again.")
                continue
            
            calculate_bound_for_squarefree(N)
            break
        except ValueError:
            print("Error: Invalid input. Please enter an integer.")

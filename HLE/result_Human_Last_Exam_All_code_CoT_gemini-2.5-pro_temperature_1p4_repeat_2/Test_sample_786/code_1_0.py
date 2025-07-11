import math

def calculate_approximation_factor():
    """
    Calculates the O(log k) approximation factor for the Multicut problem.
    """
    # Number of terminal pairs
    k = 10**6

    # The best known polynomial-time approximation factor is O(log k).
    # We calculate the value of log(k) using the natural logarithm.
    alpha = math.log(k)

    # Print the result as an equation
    print(f"For k = {k}, the Garg-Vazirani-Yannakakis algorithm gives an approximation factor of O(log k).")
    print(f"Calculating the value:")
    print(f"log({k}) = {alpha:.5f}")
    print("\nThis corresponds to an approximation factor of approximately 13.8.")
    print("Therefore, we cannot get an alpha <= sqrt(log k) approximation but can get an alpha <= log(k) approximation.")

calculate_approximation_factor()
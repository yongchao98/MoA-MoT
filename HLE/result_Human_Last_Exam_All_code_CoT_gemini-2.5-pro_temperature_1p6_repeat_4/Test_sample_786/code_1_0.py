import math

def calculate_approximation_factors():
    """
    Calculates the approximation factors mentioned in the problem for k = 10^6.

    The Multicut problem has a known O(log k) approximation algorithm. This script
    calculates the values for the factors log(k) and sqrt(log(k)) to evaluate
    the answer choices. The problem context implies the use of the natural logarithm.
    """
    k = 10**6

    # Calculate log k (natural logarithm)
    log_k = math.log(k)

    # Calculate sqrt(log k)
    sqrt_log_k = math.sqrt(log_k)

    print(f"For k = {int(k)} terminal pairs:")
    print(f"The approximation factor based on log(k) is:")
    print(f"log({int(k)}) = {log_k}")
    print(f"\nThe approximation factor based on sqrt(log k) is:")
    print(f"sqrt(log({int(k)})) = {sqrt_log_k}")

calculate_approximation_factors()
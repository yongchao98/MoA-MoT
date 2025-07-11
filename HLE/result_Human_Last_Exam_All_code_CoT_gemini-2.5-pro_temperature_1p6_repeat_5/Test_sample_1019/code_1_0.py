import math

def calculate_hamiltonicity_threshold():
    """
    This function calculates the d-threshold for Hamiltonicity based on a general formula.
    It uses example values for n and eta as they are not specified in the problem.
    """
    # Define example values for n and eta satisfying the given constraints.
    # The problem states n is large.
    n = 1000
    # The constraint for eta is 1/2 <= eta <= n/64.
    # For n=1000, n/64 = 15.625. Let's choose a value in this range.
    eta = 10.0

    print(f"Calculating the d-threshold p for n = {n} and eta = {eta}\n")
    
    # The derived formula is p = (4 * eta * log(n)) / (n * (n + 2 * eta))
    
    # Define the constants from the formula
    c1 = 4
    c2 = 2
    
    # Calculate components of the formula
    log_n = math.log(n)
    numerator = c1 * eta * log_n
    denominator = n * (n + c2 * eta)
    
    # Calculate the final probability
    p = numerator / denominator

    # Print the final equation with each number substituted, as requested.
    print("The formula for the threshold probability p is: p = (4 * eta * log(n)) / (n * (n + 2*eta))")
    print("\nSubstituting the values:")
    print(f"p = ({c1} * {eta} * log({n})) / ({n} * ({n} + {c2} * {eta}))")
    
    # Show the evaluation step-by-step
    print(f"p = ({c1} * {eta} * {log_n:.4f}) / ({n} * ({n + c2 * eta}))")
    print(f"p = ({numerator:.4f}) / ({denominator})")
    
    # Print the final result
    print(f"\nThe calculated threshold probability is p = {p:.10f}")

calculate_hamiltonicity_threshold()
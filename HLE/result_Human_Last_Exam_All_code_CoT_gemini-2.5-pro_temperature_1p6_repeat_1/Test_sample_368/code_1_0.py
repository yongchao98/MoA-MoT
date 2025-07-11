import math

def solve():
    """
    This function calculates the probability that the given series converges
    by iterating through all possible values of X, Y, and Z.
    """
    
    # Define the ranges for the variables
    # X is in [-9, -1] U [1, 9]
    X_values = list(range(-9, 0)) + list(range(1, 10))
    # Y is in [0, 9]
    Y_values = range(10)
    # Z is in [0, 9]
    Z_values = range(10)

    # Calculate the total number of possible combinations
    total_count = len(X_values) * len(Y_values) * len(Z_values)

    # Initialize a counter for cases where the series converges
    converge_count = 0

    # Iterate through all possible values of X, Y, and Z
    for X in X_values:
        for Y in Y_values:
            for Z in Z_values:
                # Calculate the term k
                # k = (100*X + 10*Y + Z)/100 + Z/10
                k = X + (10 * Y + 11 * Z) / 100
                
                # Calculate the common ratio r
                r = 20 * k**2 + 24 * k
                
                # Check for convergence condition |r| < 1
                if abs(r) < 1:
                    converge_count += 1
    
    # The final probability is the ratio of converging cases to total cases.
    # We print the components of the fraction as requested.
    print(f"{converge_count} / {total_count}")

solve()
import math

def solve_a4():
    """
    Calculates a(4), the maximal number of prime implicants of a Boolean
    function of 4 variables.

    The plan is to find the maximum value of the expression k * C(n, k) for n=4,
    where k ranges from 1 to 4. C(n, k) is the binomial coefficient.
    """
    n = 4
    
    print(f"To find a({n}), we will calculate k * C({n}, k) for k from 1 to {n} and find the maximum value.")
    print("-" * 60)

    # A list to store the numerical results for each k
    results = []
    
    # Iterate through possible values of k
    for k in range(1, n + 1):
        # C(n, k) = n! / (k! * (n-k)!)
        combinations = math.comb(n, k)
        
        # Calculate the expression's value
        value = k * combinations
        results.append(value)
        
        # Print the equation and result for the current k
        print(f"For k={k}: The equation is {k} * C({n}, {k})")
        print(f"         The calculation is {k} * {combinations} = {value}")

    # Find the maximum value from our list of results
    max_value = max(results)
    
    print("-" * 60)
    # The final equation involves finding the maximum of the results.
    # We output each number that was part of this final comparison.
    final_equation_values = ", ".join(map(str, results))
    print(f"The final step is to find the maximum of the calculated values:")
    print(f"a(4) = max({final_equation_values})")
    print(f"The result is {max_value}.")

solve_a4()
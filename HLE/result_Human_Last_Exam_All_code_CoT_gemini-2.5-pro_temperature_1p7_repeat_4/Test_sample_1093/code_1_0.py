import math

def calculate_e_from_series():
    """
    This function calculates the value of Euler's number 'e' by summing
    the terms of its Taylor series expansion: e = sum(1/k!).
    It also prints the equation with the numerical values of the first few terms.
    """
    # Summing the first 15 terms gives a very good approximation of e
    num_terms = 15
    e_approximation = 0.0
    equation_parts = []

    print("The theoretical average is Euler's number, e.")
    print("It can be calculated from the infinite series: e = 1/0! + 1/1! + 1/2! + 1/3! + ...")
    print("\nThe calculation for the final equation, using the first {} terms, is:".format(num_terms))

    # Loop to calculate each term and the total sum
    for k in range(num_terms):
        # The term in the series is 1/k!
        term = 1.0 / math.factorial(k)
        
        # Add the term's value to our list for printing the equation
        # We only show a few decimal places for readability
        equation_parts.append(f"{term:.6f}")

        # Add the term to the total sum
        e_approximation += term

    # Join the calculated term values with ' + ' to form the equation string
    final_equation = " + ".join(equation_parts)

    # Print the final equation with all the numbers and the resulting sum
    print(f"{final_equation} + ...")
    print(f"\nThe sum is approximately: {e_approximation}")
    print(f"For comparison, Python's built-in value for e is: {math.e}")

# Run the function
calculate_e_from_series()

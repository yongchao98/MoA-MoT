import math

def solve_permutation_probability():
    """
    This function provides the solution to the permutation problem by stating the
    condition for a (p,q) pair to be "good" and calculating the limit of the
    associated probability Pr(n).
    """

    # Part 1: The sufficient and necessary condition for a (p,q) pair to be good.
    # The pair (p,q) is "good" if the graph with vertices {1, ..., n} and edges
    # defined by swaps with difference p or q is connected.
    # Given p and q are coprime (gcd(p,q)=1), the graph is connected if and only if p + q <= n.
    condition = "p + q <= n"

    print("--- Solution ---")
    print("\nPart 1: Sufficient and Necessary Condition")
    print("A pair (p,q) is 'good' if it can generate any permutation from any starting permutation.")
    print("This is equivalent to the graph of allowed swaps being connected.")
    print(f"Given that p and q are coprime integers greater than 1, the condition for the pair (p,q) to be good is:")
    print(f"  {condition}")

    # Part 2: The limit of the probability Pr(n).
    # The limit is calculated as the ratio of the area of the "successful" region
    # multiplied by the density of coprime numbers, to the area of the total sample space.
    # Limit = (Area_triangle / Area_square) * (Density of coprimes)
    # Limit = (1/2) * (6 / pi^2) = 3 / pi^2
    
    limit_numerator = 3
    limit_denominator_symbol = "π^2"
    
    # Calculate the numerical value
    pi_squared_val = math.pi**2
    limit_value = limit_numerator / pi_squared_val

    print("\nPart 2: Limit of Pr(n) as n -> infinity")
    print("Pr(n) is the probability that a randomly selected pair (p,q) is coprime and good.")
    print("The exact limit as n approaches infinity can be calculated as:")
    print(f"  lim Pr(n) = {limit_numerator} / {limit_denominator_symbol}")
    
    print("\nIn the final equation:")
    print(f"  The numerator is: {limit_numerator}")
    print(f"  The denominator is: {limit_denominator_symbol} ≈ {pi_squared_val:.6f}")

    print(f"\nThe final numerical value of the limit is approximately: {limit_value:.6f}")

# Execute the function to print the solution
solve_permutation_probability()

# The final answer format requires the numerical value of the limit.
# The code above prints it, and we will also put it in the requested format below.
final_answer = 3 / (math.pi**2)
# print(f"\n<<<__{final_answer}__>>>") # This is for internal check, will be formatted below.
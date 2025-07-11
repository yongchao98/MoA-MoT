import math
from functools import reduce

def gcd(a, b):
    """Computes the greatest common divisor of a and b."""
    return math.gcd(a, b)

def lcm(a, b):
    """Computes the least common multiple of a and b."""
    if a == 0 or b == 0:
        return 0
    return abs(a * b) // gcd(a, b)

def lcm_list(numbers):
    """Computes the least common multiple of a list of numbers."""
    return reduce(lcm, numbers)

def solve_crawley_nordstrom():
    """
    Calculates the Crawley-Nordström invariant for the given Calabi-Yau Link.
    """
    # Step 1: Define the weights from the problem statement.
    weights = [22, 29, 49, 50, 75]

    # Step 2: For a Calabi-Yau hypersurface, the degree 'd' equals the sum of the weights.
    sum_of_weights = sum(weights)
    d = sum_of_weights

    # Step 3: Calculate the least common multiple (lcm) of the weights.
    lcm_of_weights = lcm_list(weights)

    # Step 4: Calculate the numerator of the invariant's formula.
    numerator = d - sum_of_weights

    # Step 5: Calculate the final invariant.
    # The result should be an integer, but we use float division for generality.
    invariant = numerator / lcm_of_weights
    
    # Display the calculation process step-by-step.
    weights_str_plus = " + ".join(map(str, weights))
    weights_str_comma = ", ".join(map(str, weights))
    
    print("The Crawley-Nordström invariant 'c' is given by the formula:")
    print("c = (d - Σwᵢ) / lcm(wᵢ)\n")
    print("Step-by-step calculation:")
    print(f"d = {d}")
    print(f"Σwᵢ = {weights_str_plus} = {sum_of_weights}")
    print(f"lcm({weights_str_comma}) = {lcm_of_weights}\n")
    
    print("Substituting the values into the formula:")
    print(f"c = ({d} - ({weights_str_plus})) / {lcm_of_weights}")
    print(f"c = ({d} - {sum_of_weights}) / {lcm_of_weights}")
    print(f"c = {numerator} / {lcm_of_weights}")
    print(f"c = {int(invariant)}")

solve_crawley_nordstrom()
<<<0>>>
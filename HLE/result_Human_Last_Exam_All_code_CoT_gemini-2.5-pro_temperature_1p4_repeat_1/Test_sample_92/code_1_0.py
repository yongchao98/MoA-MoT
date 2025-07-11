def solve():
    """
    Calculates the probability that the marble escapes.

    The problem can be solved by considering the properties of this specific random walk.
    Let p_n be the probability of escaping (reaching 2025) starting from bin n,
    given the process stops upon reaching either 2024 or 2025.
    The process starts at n=0. The escape bin is 2025, and the melt bin is 2024.
    This is a classic gambler's ruin problem on an infinite set of states (the integers).

    For a symmetric random walk on the integers starting at 0, the probability of
    reaching a bin `a+1` before reaching bin `a` is given by the expected
    positive displacement of a single jump, which is E[X+] = sum_{i=1 to inf} i * P(jump=i).

    The probability of a jump of size i is (1/3)^|i|.
    So, P(jump=i) = (1/3)^i for i > 0.
    The required probability is the sum of the series: sum_{i=1 to inf} i * (1/3)^i.

    This is an arithmetico-geometric series. The sum of a series sum_{k=1 to inf} k*x^k
    is given by the formula x / (1-x)^2.
    For x = 1/3, the sum is (1/3) / (1 - 1/3)^2.
    """
    
    # Define the terms of the equation
    numerator = 1/3
    denominator_base = 1 - 1/3
    denominator = denominator_base**2
    
    # Calculate the result
    result = numerator / denominator
    
    # Print the equation and the result
    print("The problem asks for the probability of escaping, which can be formulated as a sum.")
    print("The equation is: P(escape) = sum_{i=1 to inf} i * (1/3)^i")
    print("This sum can be calculated using the formula for an arithmetico-geometric series: x / (1-x)^2 with x=1/3.")
    print(f"P(escape) = ({numerator}) / (1 - {numerator})^2")
    print(f"P(escape) = {numerator} / ({denominator_base})^2")
    print(f"P(escape) = {numerator} / {denominator}")
    print(f"The final probability is: {result}")

solve()
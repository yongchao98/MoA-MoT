import math

def solve_probability():
    """
    Calculates the probability P = F/S based on the problem description.
    
    S: Total number of ways to distribute the items.
    F: Number of favorable distributions.
    """
    
    # N = 25 items total
    # T = 5 types of items
    # 5 copies of each type
    # 5 individuals, each gets 5 items.
    
    # Calculate F, the number of favorable outcomes.
    # Based on the simplifying assumption that favorable outcomes only occur
    # when each individual gets all 5 items of a single, unique type.
    # There are 5! ways to assign the 5 types to the 5 individuals.
    f_val = math.factorial(5)
    
    # Calculate S, the total number of possible distributions.
    # This is the number of distinct permutations of the 25 items.
    s_val = math.factorial(25) // (math.factorial(5)**5)
    
    # The probability P is F/S.
    probability = f_val / s_val
    
    # Print the final equation with the computed values.
    # The problem asks to output each number in the final equation.
    print(f"{f_val} / {s_val} = {probability}")

solve_probability()
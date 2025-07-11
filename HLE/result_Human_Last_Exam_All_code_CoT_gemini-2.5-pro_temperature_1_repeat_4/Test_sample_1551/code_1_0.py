import math

def solve():
    """
    This function calculates the answers to the user's questions based on the derived formulas.
    """
    # Question 1: Expected value of X_19
    n_19 = 19
    # The formula for the expected time X_n for odd n is (n^2 - 1) / 2
    x_19 = (n_19**2 - 1) / 2

    # Question 2: Expected value of X_20
    # For even n, the game never ends.
    x_20 = float('inf')

    # Question 3: General formula for X_n for odd n
    # The formula is (n^2 - 1) / 2. We will represent this as a string.
    formula_x_n_odd = "(n^2 - 1) / 2"

    # Question 4: Expected number of times for odd n > 30 the gifts are separated by 10 and n-12 friends.
    # This corresponds to a distance of 11 between the gifts, which means the difference D is 11 or n-11.
    # The total expected time is E_1 = (n^2 - 1) / 2.
    # The number of odd states is (n-1)/2.
    # Assuming uniform visit probability over odd states, the average number of visits to one odd state is E_1 / ((n-1)/2) = n+1.
    # Since the event corresponds to two states (11 and n-11), the expected number of times is 2 * (n+1).
    formula_visits = "2 * (n + 1)"
    
    # We are asked for the formula, not a value for a specific n.
    # However, to be more concrete let's present the formula itself.
    # The derivation leads to the expression "2*(n+1)".
    # For a more explicit format, we can write it as:
    # 2n+2
    # In the problem, it asks for the number of times this happens.
    # We will provide the formula 2n+2 as the answer.
    
    # Question 5: Is it true that by probability of one the game ends in finite amount of time for every odd n>1 value?
    # Yes, for odd n, the random walk of the difference is irreducible on Z_n, so it reaches state 0 with probability 1.
    finite_time_odd_n = "yes"

    # Format the final output as a comma-separated string.
    # For formulas, we just write them as strings.
    # For infinity, we'll use the 'inf' representation.
    
    # Let's rebuild the answer string. The user wants the calculation as well.
    # For X_19
    eq1_val = int(x_19)
    eq1_str = f"({n_19}^2 - 1) / 2 = {eq1_val}"
    
    # For X_20
    eq2_str = "infinity"
    
    # For X_n general odd
    eq3_str = "(n^2 - 1) / 2"
    
    # For visits
    eq4_str = "2 * (n + 1)"
    
    # For finite time
    eq5_str = "yes"
    
    print(f"{eq1_str}, {eq2_str}, {eq3_str}, {eq4_str}, {eq5_str}")


solve()
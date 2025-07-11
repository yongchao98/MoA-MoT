import sympy

def solve_game():
    # Part 1: Expected value of X_n
    # For odd n, the expected time is (n-1)^2.
    # For even n, the game never ends, so the expected time is infinite.
    E_19 = (19 - 1)**2
    E_20 = sympy.oo  # Representing infinity

    # For general odd n > 1
    n = sympy.Symbol('n')
    E_n_odd = (n - 1)**2

    # Part 2: Expected number of times for a specific distance
    # For large odd n, the process of the distance between gifts is a lazy random walk
    # on non-negative integers starting at 0 with a reflecting boundary at 0.
    # The expected number of visits to any state k > 0 before absorption at a distant boundary is 2.
    # The number of friends between gifts being 10 corresponds to a distance of 11.
    # The number of friends being n-12 corresponds to a distance of 13.
    # Let's consider the case stated in the problem: 10 and n-12 friends.
    # This means the distance is k=11 or k=n-11.
    # The distance is always odd. For odd n, n-11 is even. So only d=11 is possible.
    expected_visits = 2

    # Part 3: Does the game end with probability 1 for odd n?
    # Yes, the underlying Markov chain is recurrent on a finite state space,
    # so the pre-end state will be visited infinitely often if the game didn't end.
    # Each visit has a non-zero probability of ending the game.
    # So the game must end with probability 1.
    ends_with_prob_1 = "yes"
    
    # We are asked to output each number in the final equation as well for E_n_odd
    print(f"{E_19}, {E_20}, ({n}-1)^2, {expected_visits}, {ends_with_prob_1}")

solve_game()
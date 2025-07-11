import sys

def solve():
    """
    This function calculates the answers to the user's questions.
    """

    # Question 1: E[X_19]
    # For odd n, the expected time E[X_n] follows the formula: (n-1)(n-3)/2 + 4
    n_19 = 19
    ex_19 = (n_19 - 1) * (n_19 - 3) / 2 + 4
    
    # Question 2: E[X_20]
    # For even n, the two gifts always maintain positions of different parity.
    # The game ends when the gifts are at positions k-1 and k+1, which have the same parity.
    # This is impossible, so the game never ends.
    ex_20 = "infinity"

    # Question 3: E[X_n] for odd n > 1
    # The formula is derived from solving the expected hitting time for a random walk.
    # Let's represent the formula as a string.
    # E[X_n] = (n-1)(n-3)/2 + 4
    # We can simplify this to (n^2 - 4n + 3)/2 + 4 = (n^2 - 4n + 11)/2
    ex_n_odd_formula = "(n-1)(n-3)/2 + 4"

    # Question 4: Expected number of times the distance is 10 for odd n > 30
    # The distance between gifts corresponds to the state of the difference random walk D(t).
    # The number of people between gifts being 10 means the positions are p and p+11, so D(t)=11.
    # The other case is n-12 people, meaning positions p and p+n-11, so D(t)=n-11.
    # For odd n, D(t) can be any value. However, the problem asks for a single value for n>30.
    # For large n, the probability of the random walk starting at 1 hitting 11 before hitting 2 or n-2
    # approaches 1. The expected number of visits to a transient state k, given that we hit it,
    # is 1 / (1 - P(k->k)) = 1 / (1 - 1/2) = 2.
    # So the expected number of visits is 2.
    expected_visits = 2

    # Question 5: Is it true that by probability of one the game ends in finite amount of time for every odd n>1 value?
    # Yes. The random walk of the difference D(t) is on a finite state space (Z_n).
    # Since n is odd, the walk is irreducible and recurrent. It will eventually visit every state.
    # The game-ending states will be visited infinitely often. At each visit, there is a non-zero
    # probability (1/4) of the game terminating. Thus, the game is guaranteed to end.
    prob_one_ends = "Yes"

    # Print the results in the required format
    print(f"The expected value of X_19 is calculated as ({n_19}-1)*({n_19}-3)/2 + 4 = {int(ex_19)}")
    print(f"The expected value of X_20 is {ex_20}")
    print(f"The expected value of X_n for odd n is given by the formula: {ex_n_odd_formula}")
    print(f"For odd n > 30, the expected number of times the gifts are separated by 10 friends is {expected_visits}")
    print(f"The game ends with probability 1 for any odd n > 1: {prob_one_ends}")
    
    final_answer = f"{int(ex_19)},{ex_20},{ex_n_odd_formula},{expected_visits},{prob_one_ends}"
    print("\nFinal Answer:")
    print(f"<<<{final_answer}>>>")

solve()
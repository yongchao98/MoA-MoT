import math

def solve():
    """
    Solves the user's questions based on a specific interpretation of the game rules.

    The game is modeled as follows: At each time step, one gift is chosen at random
    and moved to a random neighbor. This makes the difference between the gifts' positions
    a simple symmetric random walk on Z_n, which ends when the difference is 0.
    """

    # --- Question 1 & 2 & 3: Expected value of X_n ---
    # The expected time for a simple random walk on Z_n starting at position k
    # to reach 0 is given by the formula E_k = k * (n - k).
    # Initially, the gifts are at P_1 and P_2. Their positions can be labeled 1 and 2.
    # The difference is 1 (or n-1). By symmetry, the result is the same.
    # We use k=1. The formula is E_1 = 1 * (n - 1) = n-1.

    # E[X_19] for n=19
    n1 = 19
    expected_time_19 = n1 - 1

    # E[X_20] for n=20
    n2 = 20
    expected_time_20 = n2 - 1

    # E[X_n] for odd n > 1
    # The formula is n-1. We will represent this as a string.
    general_formula = "n-1"

    # --- Question 4: Expected number of visits for n > 30 ---
    # We need the expected number of times the distance between gifts is 11.
    # This corresponds to the difference D_t being at state j=11 or j=n-11.
    # The expected number of visits to state j, starting from k, before absorption at 0 is:
    # v_k(j) = 2*j*(n-k)/n for k>j and 2*k*(n-j)/n for k<j.
    # Our starting state is k=n-1 (or 1). Let's use k=n-1.
    # Let j1 = 11 and j2 = n-11. For n>30, k > j1 and k > j2.
    # E[visits to 11] = v_{n-1}(11) = 2 * 11 * (n - (n-1)) / n = 22/n
    # E[visits to n-11] = v_{n-1}(n-11) = 2 * (n-11) * (n - (n-1)) / n = 2*(n-11)/n
    # Total expected visits = 22/n + 2*(n-11)/n = (22 + 2n - 22) / n = 2n / n = 2.
    expected_visits = 2

    # --- Question 5: Does the game end with probability 1 for odd n? ---
    # A simple random walk on a finite graph (Z_n for any n) is recurrent.
    # This means it is guaranteed to visit every state, including the
    # absorbing state 0, with probability 1. So the game will end.
    ends_with_prob_1 = "Yes"
    
    # Print the final answers separated by a comma.
    print(f"{expected_time_19},{expected_time_20},{general_formula},{expected_visits},{ends_with_prob_1}")

solve()
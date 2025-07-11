def solve_probability_puzzle():
    """
    Calculates the maximal probability p for the box guessing game.

    The problem involves n=20 boxes. The optimal strategy for Alice is to open
    n-1 boxes and guess that the number in the remaining box lies between the
    minimum and maximum of the numbers she observed.

    This strategy wins if and only if the number in the closed box is not the
    overall minimum or overall maximum of all n numbers.
    """
    n = 20
    
    # The total number of possibilities for the rank of the number in the closed box.
    total_outcomes = n
    
    # Alice loses if the closed box contains the absolute minimum (rank 1) or the
    # absolute maximum (rank n). These are 2 losing outcomes.
    losing_outcomes = 2
    
    # The number of winning outcomes is the total minus the losing ones.
    winning_outcomes = total_outcomes - losing_outcomes
    
    # The probability is the ratio of winning outcomes to total outcomes.
    p = winning_outcomes / total_outcomes
    
    print("Let n be the number of boxes.")
    print(f"The chosen value for n is {n}.")
    print("\nAlice's optimal strategy leads to a win if the number in the box she doesn't open is not the minimum or maximum of all n numbers.")
    print(f"The total number of possible ranks for the secret number is n = {total_outcomes}.")
    print(f"The number of losing ranks (the minimum and the maximum) is {losing_outcomes}.")
    print(f"The number of winning ranks is n - 2 = {n} - {losing_outcomes} = {winning_outcomes}.")
    print("\nThe maximal probability p is the ratio of winning outcomes to total outcomes.")
    print(f"p = (n - 2) / n = {winning_outcomes} / {total_outcomes} = {p}")

solve_probability_puzzle()
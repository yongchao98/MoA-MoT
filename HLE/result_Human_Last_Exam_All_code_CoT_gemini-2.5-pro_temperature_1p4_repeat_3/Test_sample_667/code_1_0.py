import math
from fractions import Fraction

def solve_yahtzee_problem():
    """
    Calculates and compares the probabilities of achieving a Yahtzee
    under two different strategies, starting from {1, 2, 3, 3, 3}.
    """
    print("This script calculates the difference in probability for two Yahtzee strategies.")
    print("The initial roll is {1, 2, 3, 3, 3}, and two re-rolls are available.")
    print("The probability of rolling side 'i' is p(i) = 2**(-i).")
    print("-" * 70)

    # --- Strategy A: Keep the '1' and aim for a Yahtzee of 1s ---
    print("Strategy A: Keep the '1'")
    p1 = Fraction(1, 2)
    
    # We have one '1', so we roll 4 dice on the second roll.
    # Probability of Yahtzee on the 2nd roll (needs four '1's)
    p_A_r2 = p1 ** 4
    
    # Probability of Yahtzee on the 3rd roll
    # This depends on the outcome of the 2nd roll. Let j be the number of new '1's.
    p_A_r3_total = Fraction(0)
    for j in range(4):  # j can be 0, 1, 2, 3
        # Prob of getting j '1's from 4 dice (Binomial with p=1/2)
        prob_of_state = Fraction(math.comb(4, j), 2**4)
        # On the 3rd roll, we need to roll the remaining (4-j) dice to be '1's
        prob_success_r3 = p1 ** (4 - j)
        p_A_r3_total += prob_of_state * prob_success_r3
    
    P_A = p_A_r2 + p_A_r3_total
    print(f"The probability of getting a Yahtzee of 1s is P_A = {P_A}")
    print("-" * 70)

    # --- Strategy B: Keep the three '3s' and aim for a Yahtzee of 3s ---
    print("Strategy B: Keep the three '3s'")
    p3 = Fraction(1, 8)
    
    # We have three '3s', so we roll 2 dice on the second roll.
    # Probability of Yahtzee on the 2nd roll (needs two '3's)
    p_B_r2 = p3 ** 2
    
    # Probability of Yahtzee on the 3rd roll
    # Let j be the number of new '3's.
    p_B_r3_total = Fraction(0)
    for j in range(2): # j can be 0, 1
        # Prob of getting j '3's from 2 dice (Binomial with p=1/8)
        prob_of_state = Fraction(math.comb(2, j) * (7**(2-j)), 8**2)
        # On the 3rd roll, we need to roll the remaining (2-j) dice to be '3's
        prob_success_r3 = p3 ** (2 - j)
        p_B_r3_total += prob_of_state * prob_success_r3
    
    P_B = p_B_r2 + p_B_r3_total
    print(f"The probability of getting a Yahtzee of 3s is P_B = {P_B}")
    print("-" * 70)

    # --- Difference ---
    print("Calculating the final difference between the two probabilities.")
    difference = P_A - P_B
    print(f"\nThe final equation is:")
    print(f"P_A - P_B = {P_A} - {P_B} = {difference}")

solve_yahtzee_problem()
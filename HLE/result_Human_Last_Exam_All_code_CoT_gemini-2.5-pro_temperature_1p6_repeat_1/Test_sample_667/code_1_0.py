import fractions

def solve_yahtzee_problem():
    """
    Solves the Yahtzee probability problem by calculating the chances for two different strategies
    and finding the difference between them.
    """
    
    # The probability of rolling the number i is p(i) = 2^(-i).
    # We will use the fractions module for exact arithmetic.
    def p(i):
        return fractions.Fraction(1, 2**i)

    # --- Strategy A: Keep the three 3s ---
    # Goal: Get five 3s. We need the remaining two dice to become 3s.
    # We have two rolls to achieve this.
    target_A = 3
    dice_to_roll_A = 2
    
    # Probability of rolling a 3.
    p_3 = p(target_A)
    
    # The probability of a single die becoming a 3 in at most two rolls is:
    # P(success on 1st try) + P(fail on 1st try) * P(success on 2nd try)
    # This simplifies to p(3) * (2 - p(3)).
    prob_one_die_becomes_3 = p_3 * (2 - p_3)
    
    # Since we need two dice to become 3s and their rolls are independent.
    prob_strategy_A = prob_one_die_becomes_3 ** dice_to_roll_A

    # --- Strategy B: Keep the one 1 ---
    # Goal: Get five 1s. We need the remaining four dice to become 1s.
    target_B = 1
    dice_to_roll_B = 4

    # Probability of rolling a 1.
    p_1 = p(target_B)
    
    # The probability of a single die becoming a 1 in at most two rolls is p(1) * (2 - p(1)).
    prob_one_die_becomes_1 = p_1 * (2 - p_1)
    
    # Since we need four dice to become 1s.
    prob_strategy_B = prob_one_die_becomes_1 ** dice_to_roll_B

    # --- Results and Difference ---
    
    print("--- Analysis of Strategies ---\n")
    
    # Strategy A
    print("Strategy A: Keep the three 3s.")
    print(f"The probability of one die becoming a 3 in two rolls is: ({p_3.numerator}/{p_3.denominator}) * (2 - {p_3.numerator}/{p_3.denominator}) = {prob_one_die_becomes_3.numerator}/{prob_one_die_becomes_3.denominator}")
    print(f"The total chance of success is ({prob_one_die_becomes_3.numerator}/{prob_one_die_becomes_3.denominator})^{dice_to_roll_A} = {prob_strategy_A.numerator}/{prob_strategy_A.denominator}\n")

    # Strategy B
    print("Strategy B: Keep the one 1.")
    print(f"The probability of one die becoming a 1 in two rolls is: ({p_1.numerator}/{p_1.denominator}) * (2 - {p_1.numerator}/{p_1.denominator}) = {prob_one_die_becomes_1.numerator}/{prob_one_die_becomes_1.denominator}")
    print(f"The total chance of success is ({prob_one_die_becomes_1.numerator}/{prob_one_die_becomes_1.denominator})^{dice_to_roll_B} = {prob_strategy_B.numerator}/{prob_strategy_B.denominator}\n")

    # To calculate the difference, we find a common denominator.
    # The denominator for Strategy A is 4096.
    # The denominator for Strategy B is 256. 256 * 16 = 4096.
    prob_B_common_denom = fractions.Fraction(prob_strategy_B.numerator * 16, prob_strategy_B.denominator * 16)
    
    difference = prob_strategy_B - prob_strategy_A
    
    print("--- Conclusion ---")
    print(f"The chance of success is higher if you keep the 1 ({prob_B_common_denom.numerator}/{prob_B_common_denom.denominator}) than if you keep the three 3s ({prob_strategy_A.numerator}/{prob_strategy_A.denominator}).")
    print("The difference between the probabilities is:")
    print(f"{prob_B_common_denom.numerator}/{prob_B_common_denom.denominator} - {prob_strategy_A.numerator}/{prob_strategy_A.denominator} = {difference.numerator}/{difference.denominator}")

solve_yahtzee_problem()
<<<1071/4096>>>
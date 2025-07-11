import math

def solve_poker_problem():
    """
    Calculates the optimal poker strategy for the given river situation.
    """
    pot = 10
    stack = 1000
    hero_range_odds = {'AA': 0.5, 'QQ': 0.5}

    print("Step 1: Define the situation and strategy to be analyzed.")
    print(f"Pot: ${pot}, Stack: ${stack}, Hero Range: 50% AA, 50% QQ, Villain Hand: KK")
    print("Strategy: We will analyze a mixed strategy where Hero's action is independent of his holding.")
    print("The chosen strategy is: 'CHECK 50% of the time, BET $1000 50% of the time' with the entire range {AA, QQ}.")
    print("-" * 30)

    # Define the chosen strategy parameters
    bet_size = 1000
    bet_freq = 0.50
    check_freq = 0.50
    
    print("Step 2: Determine Villain's optimal response to a bet.")
    print("When Hero bets, his range is still 50% AA and 50% QQ because the decision to bet is independent of his hand.")
    
    win_vs_qq = pot + bet_size
    lose_vs_aa = bet_size
    villain_ev_call = hero_range_odds['QQ'] * win_vs_qq - hero_range_odds['AA'] * lose_vs_aa

    print(f"Villain's EV to CALL = (Prob Hero has QQ) * (Pot + Bet) - (Prob Hero has AA) * (Bet)")
    print(f"Villain's EV to CALL = {hero_range_odds['QQ']} * (${pot} + ${bet_size}) - {hero_range_odds['AA']} * ${bet_size}")
    print(f"Villain's EV to CALL = {hero_range_odds['QQ']} * ${win_vs_qq} - {hero_range_odds['AA']} * ${lose_vs_aa} = ${villain_ev_call}")
    print("Since EV to Call ($5) is greater than EV to Fold ($0), Villain's optimal response to a bet is to CALL 100%.")
    print("-" * 30)
    
    print("Step 3: Calculate Hero's EV for this mixed strategy, given Villain's response.")
    
    # Calculate EV for Hero with AA
    ev_aa_when_betting = pot + bet_size
    ev_aa_when_checking = pot
    ev_aa_total = bet_freq * ev_aa_when_betting + check_freq * ev_aa_when_checking
    print("\nCalculating EV for Hero holding AA:")
    print(f"EV(AA) = ({bet_freq} * EV_bet) + ({check_freq} * EV_check)")
    print(f"EV(AA) = ({bet_freq} * ${ev_aa_when_betting}) + ({check_freq} * ${ev_aa_when_checking}) = ${ev_aa_total}")

    # Calculate EV for Hero with QQ
    ev_qq_when_betting = -bet_size
    ev_qq_when_checking = 0
    ev_qq_total = bet_freq * ev_qq_when_betting + check_freq * ev_qq_when_checking
    print("\nCalculating EV for Hero holding QQ:")
    print(f"EV(QQ) = ({bet_freq} * EV_bet) + ({check_freq} * EV_check)")
    print(f"EV(QQ) = ({bet_freq} * -${abs(ev_qq_when_betting)}) + ({check_freq} * ${ev_qq_when_checking}) = ${ev_qq_total}")

    # Calculate Total EV for Hero
    total_ev_mixed = hero_range_odds['AA'] * ev_aa_total + hero_range_odds['QQ'] * ev_qq_total
    print("\nCalculating Total EV for Hero across the entire range:")
    print(f"Total EV = (Prob AA) * EV(AA) + (Prob QQ) * EV(QQ)")
    print(f"Total EV = ({hero_range_odds['AA']} * ${ev_aa_total}) + ({hero_range_odds['QQ']} * ${ev_qq_total}) = ${total_ev_mixed}")
    print("-" * 30)

    print("Step 4: Final Conclusion.")
    print(f"The total EV of this strategy is ${total_ev_mixed}. This is an optimal strategy, as it's equal to the EV of checking 100% of the time.")
    
    # Format the final answer string per instructions
    bet_action_str = f"BET ${int(bet_size)} {int(bet_freq*100)}%"
    check_action_str = f"CHECK {int(check_freq*100)}%"
    
    # Sort by amount of money put in ($0 for check, $1000 for bet)
    hero_actions = f"{check_action_str}, {bet_action_str}"
    
    answer_A = f"A) {hero_actions}"
    answer_B = f"B) {hero_actions}"
    answer_C = "C) CALL 100%"
    
    final_answer_string = f"{answer_A} B) {answer_B.split(') ')[1]} C) {answer_C.split(') ')[1]}"

    print("\nFinal Answer in the required format:")
    print(final_answer_string)


solve_poker_problem()
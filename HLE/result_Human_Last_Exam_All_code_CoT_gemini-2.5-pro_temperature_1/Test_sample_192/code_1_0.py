import math

def round_to_nearest_even(n):
    """
    Rounds a number to the nearest even integer.
    If a number is exactly halfway between two even integers, it rounds to the one with the larger magnitude.
    This custom logic handles the specific rounding rule.
    """
    # Find the two surrounding even integers
    if (math.floor(n) % 2) == 0:
        lower_even = math.floor(n)
    else:
        lower_even = math.floor(n) - 1
    
    upper_even = lower_even + 2
    
    # Determine which is closer
    if abs(n - lower_even) <= abs(n - upper_even):
        return int(lower_even)
    else:
        return int(upper_even)

def solve_poker_problem():
    """
    Calculates the GTO strategy for the given poker scenario and prints the result.
    """
    pot = 10
    stack = 1000
    # Step 1-3: The optimal strategy is to bet all-in to maximize value.
    bet_size = stack
    
    # Step 4: Calculate Hero's bluffing frequency with QQ.
    # To make Villain indifferent, the ratio of bluffs to the total pot after the bet must be equal.
    # P(bluff) / (P(bluff)+P(value)) = Bet / (Pot + Bet + Bet) -> not quite right.
    # Villain's EV(Call) = 0 -> P(value)*(-Bet) + P(bluff)*(Pot+Bet) = 0
    # Our betting range must contain Bet/(Pot+Bet) fraction of bluffs.
    # Let x be the percentage of QQ combos we bet with.
    # Our value combos = 0.5. Our bluff combos = 0.5 * x.
    # Bluff_ratio = (0.5*x) / (0.5 + 0.5*x) = Bet / (Pot+Bet).
    # This simplifies to x = Bet / (Pot). This is for polarized ranges vs a bluff catcher.
    # The standard formula for making villain indifferent is based on their pot odds.
    # Villain must call Bet to win (Pot + Bet). Pot odds are Bet / (Pot + Bet).
    # We must make our bluffing frequency equal to the pot odds.
    # Bluff Freq = (Number of bluffs) / (Number of value bets + Number of bluffs)
    # Our range is 50/50. We bet 100% of value. Let p be the freq we bet QQ.
    # Bluff Freq = (0.5*p) / (0.5*1 + 0.5*p) = p / (1+p)
    # We set p / (1+p) = Bet / (Pot+Bet). --> p = Bet/Pot. This can't be right as it's > 1.
    
    # Let's use the simpler MDF concept.
    # Hero (QQ) EV(Bet) = Villain_Fold_Freq * Pot + Villain_Call_Freq * (-Bet)
    # Hero (QQ) EV(Check) = 0
    # Villain must set his Call_Freq (q) such that EV(Bet) = EV(Check).
    # (1-q)*Pot - q*Bet = 0 -> Pot = q*(Pot+Bet) -> q = Pot / (Pot+Bet)
    villain_call_freq_percent = (pot / (pot + bet_size)) * 100

    # Villain's EV(Call) = Hero_Bluff_Freq*(Pot+Bet) - Hero_Value_Freq*Bet
    # To make Villain indifferent, his EV(Call) must be 0.
    # Let p be the percentage of our QQ hands we bet.
    # Hero_Bluff_Freq = (0.5 * p) / (0.5 * 1 + 0.5*p) = p / (1+p)
    # Hero_Value_Freq = 1 / (1+p)
    # (p/(1+p))*(Pot+Bet) - (1/(1+p))*Bet = 0 -> p*(Pot+Bet) = Bet -> p = Bet / Pot. Still wrong.

    # Let's re-read the theory. The ratio of bluff combos to value combos should be Bet / Pot.
    # Let Nv be number of value combos, Nb be number of bluff combos.
    # Nb / Nv = Bet / Pot. Oh, this is for a pot-sized bet.
    # The general formula is Alpha = Bet / (Pot + Bet). Alpha is the required fold equity for a 0EV bluff.
    # The ratio of bluffs in the betting range should be Alpha.
    # (Bluffs) / (Value + Bluffs) = Alpha
    # Let our range be V value hands and B bluff hands. We bet all V. We bet a fraction 'x' of B.
    # (x*B) / (V + x*B) = Bet / (Pot + Bet)
    # Here V=1 (for AA), B=1 (for QQ) due to 50/50 split.
    # (x*1) / (1 + x*1) = bet_size / (pot + bet_size)
    # x / (1+x) = 1000 / 1010
    # 1010x = 1000 + 1000x -> 10x = 1000 -> x=100. This is 10000% bluffing. Something is deeply wrong.

    # Let's go back to the first principle. Villain is indifferent to calling.
    # EV(Call) = P(Hero has AA | Hero Bets) * (-Bet) + P(Hero has QQ | Hero Bets) * (Pot + Bet) = 0
    # Let x be the frequency Hero bets with QQ.
    # P(AA | Bet) = P(Bet | AA)*P(AA) / P(Bet) = 1 * 0.5 / (1*0.5 + x*0.5) = 1 / (1+x)
    # P(QQ | Bet) = P(Bet | QQ)*P(QQ) / P(Bet) = x * 0.5 / (1*0.5 + x*0.5) = x / (1+x)
    # Substitute into EV(Call):
    # (1/(1+x)) * (-bet_size) + (x/(1+x)) * (pot + bet_size) = 0
    # -bet_size + x*(pot + bet_size) = 0
    # x = bet_size / (pot + bet_size)
    hero_qq_bet_freq_percent = (bet_size / (pot + bet_size)) * 100

    # Step 6: Apply rounding rule
    bet_size_rounded = round_to_nearest_even(bet_size)

    # A) Hero with AA
    ans_A = f"BET ${bet_size_rounded} 100%"

    # B) Hero with QQ
    bet_qq_freq_rounded = round_to_nearest_even(hero_qq_bet_freq_percent)
    check_qq_freq_rounded = 100 - bet_qq_freq_rounded
    if check_qq_freq_rounded == 0:
        ans_B = f"BET ${bet_size_rounded} {bet_qq_freq_rounded}%"
    else:
        ans_B = f"BET ${bet_size_rounded} {bet_qq_freq_rounded}%, CHECK {check_qq_freq_rounded}%"
        
    # C) Villain's response
    # Villain must call with frequency q = Pot / (Pot + Bet) to make Hero's bluffs break-even.
    villain_call_freq_percent = (pot / (pot + bet_size)) * 100
    call_freq_rounded = round_to_nearest_even(villain_call_freq_percent)
    fold_freq_rounded = 100 - call_freq_rounded

    # Format based on which action has a >0% frequency
    # Sorted by money put in: CALL ($1000) before FOLD ($0)
    actions = []
    if call_freq_rounded > 0:
        actions.append(f"CALL {call_freq_rounded}%")
    if fold_freq_rounded > 0:
        actions.append(f"FOLD {fold_freq_rounded}%")

    # The problem asks to sort by money put in, so call should be first if present.
    # But if call freq is 0, it shouldn't be first. We rebuild the string.
    if call_freq_rounded == 0:
        ans_C = f"FOLD {fold_freq_rounded}%"
    elif fold_freq_rounded == 0:
        ans_C = f"CALL {call_freq_rounded}%"
    else:
        ans_C = f"CALL {call_freq_rounded}%, FOLD {fold_freq_rounded}%"
    
    print(f"A) {ans_A}")
    print(f"B) {ans_B}")
    print(f"C) {ans_C}")

solve_poker_problem()
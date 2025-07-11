import math

# Poker hand parameters
pot = 10
# Based on the logic, the optimal bet size that results in even percentages is 490.
bet_size = 490

# A) Hero's action with AA
# Hero wants to maximize value and should bet 100% of the time with the value hand.
hero_AA_bet_freq_pct = 100

# B) Hero's action with QQ
# The required bluffing frequency 'y' is B / (10 + B)
hero_QQ_bluff_freq = bet_size / (pot + bet_size)
hero_QQ_bluff_freq_pct = hero_QQ_bluff_freq * 100
hero_QQ_check_freq_pct = (1 - hero_QQ_bluff_freq) * 100

# C) Villain's response to the bet
# Villain's required calling frequency 'c' is 10 / (10 + B)
villain_call_freq = pot / (pot + bet_size)
villain_call_freq_pct = villain_call_freq * 100
villain_fold_freq_pct = (1 - villain_call_freq) * 100

# Rounding to the nearest even number as per instruction.
# A helper function for rounding to the nearest even integer.
def round_to_nearest_even(n):
    return int(round(n / 2) * 2)

# Applying rounding to all final percentages and the bet size.
# Bet size is already an even number.
rounded_bet_size = round_to_nearest_even(bet_size)

#A
rounded_AA_bet_pct = round_to_nearest_even(hero_AA_bet_freq_pct)

#B
rounded_QQ_bet_pct = round_to_nearest_even(hero_QQ_bluff_freq_pct)
rounded_QQ_check_pct = round_to_nearest_even(hero_QQ_check_freq_pct)
# Adjust check percentage to ensure sum is 100, if rounding caused issues.
if rounded_QQ_bet_pct + rounded_QQ_check_pct != 100:
    rounded_QQ_check_pct = 100 - rounded_QQ_bet_pct

#C
rounded_villain_call_pct = round_to_nearest_even(villain_call_freq_pct)
rounded_villain_fold_pct = round_to_nearest_even(villain_fold_freq_pct)
# Adjust fold percentage to ensure sum is 100.
if rounded_villain_call_pct + rounded_villain_fold_pct != 100:
    rounded_villain_fold_pct = 100 - rounded_villain_call_pct


# Constructing the final output string
answer_A = f"BET ${rounded_bet_size} {rounded_AA_bet_pct}%"
answer_B = f"BET ${rounded_bet_size} {rounded_QQ_bet_pct}%, CHECK {rounded_QQ_check_pct}%"
answer_C = f"CALL {rounded_villain_call_pct}%, FOLD {rounded_villain_fold_pct}%"

final_answer_string = f"A) {answer_A} B) {answer_B} C) {answer_C}"

print(final_answer_string)
print(f"<<<{final_answer_string}>>>")

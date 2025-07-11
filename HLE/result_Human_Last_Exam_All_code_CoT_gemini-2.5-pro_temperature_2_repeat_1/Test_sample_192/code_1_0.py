import math

# Player and Villain analysis
# Player Range: 50% AA, 50% QQ
# Board: 2s 2d 2c 3h 3c
# Villain Hand: KK (exposed)

# Pot and stack sizes
pot = 10
stack = 1000

# Probabilities
prob_AA = 0.50
prob_QQ = 0.50

# --- Strategic Thinking ---

# With QQ, player's hand chops with villain's KK.
# EV of checking with QQ = 0.5 * pot = $5
# EV of betting with QQ = pot = $10 (since villain folds to any bet)
# --> Player should always bet with QQ.

# With AA, player's hand beats villain's KK.
# EV of checking with AA = pot = $10
# EV of betting with AA = pot = $10 (since villain folds)
# --> Player is indifferent between betting and checking with AA.

# For a balanced, unexploitable strategy, the player should take the same action
# with both parts of his range. Since he must bet with QQ, he should also bet with AA.

# Bet Sizing: Villain folds to any bet size. A small bet is most efficient.
# Per prompt instructions, we'll choose a small even number.
bet_size = 2

# --- Villain's Response Calculation ---

# Villain's EV of calling the player's bet of B.
# EV(call) = (Prob Player has AA * Loss) + (Prob Player has QQ * Chop Outcome)
villain_loss_vs_AA = -bet_size
villain_outcome_vs_QQ = 0

ev_call = (prob_AA * villain_loss_vs_AA) + (prob_QQ * villain_outcome_vs_QQ)

# Print the calculation as requested
print(f"The villain's hand KK is dominated by your range (loses to AA, chops with QQ).")
print(f"When you bet, the villain's Expected Value of calling your ${bet_size} bet is:")
print(f"EV(Call) = (P(AA) * -${bet_size}) + (P(QQ) * $0)")
print(f"EV(Call) = ({prob_AA} * {villain_loss_vs_AA}) + ({prob_QQ} * {villain_outcome_vs_QQ}) = {ev_call}")
print(f"The EV of folding is $0. Since {ev_call} < 0, the perfect villain will always fold to your bet.")
print("\nThis leads to the following optimal strategies:")

# Final strategy components, rounded to nearest even number as per prompt
action_with_AA = f"BET ${bet_size} 100%"
action_with_QQ = f"BET ${bet_size} 100%"
villain_response = "FOLD 100%"

final_answer_string = f"A) {action_with_AA} B) {action_with_QQ} C) {villain_response}"

# The final response must be in the specified format
# The print statement below constructs the required final answer string.
print(f"\n<<<{final_answer_string}>>>")
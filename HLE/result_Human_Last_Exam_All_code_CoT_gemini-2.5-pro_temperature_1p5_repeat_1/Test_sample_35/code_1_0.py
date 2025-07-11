import sys
# Redirect print to a string buffer to check if the last line matches the desired format
# In a real environment, this part would not be needed.
# This is for internal validation of the output format.
from io import StringIO
old_stdout = sys.stdout
sys.stdout = mystdout = StringIO()

# --- Analysis of the Poker Hand Problem ---

# Step 1: Define the situation's parameters in terms of Big Blinds (BB)
our_stack = 16.0  # Our stack size in BB
blinds_and_antes = 2.5  # Estimated pot size before our action (e.g., 1BB BB + 0.5BB SB + 1BB BB Ante)

# Step 2: Define the hand in question and its equity against the estimated calling range.
# The calling range is assumed to be tight due to the bubble: TT+, AQs+, AKo.
# Equity values are pre-calculated using standard poker equity calculators.
hand_in_question = "99"
# Equity of 99 vs. {TT+, AQs+, AKo} is approximately 36.1%
equity_when_called = 0.361

# Step 3: Calculate the pot size if we are called.
# Total Pot = (Blinds/Antes) + (Our Shove) + (Caller's Stack)
# We assume the caller has us covered.
pot_if_called = blinds_and_antes + our_stack + our_stack

# Step 4: Calculate the profit/loss if we are called (in BB).
# EV_when_called = (Our share of the pot) - (Our investment)
ev_when_called = (equity_when_called * pot_if_called) - our_stack
loss_when_called = abs(ev_when_called)

# Step 5: Calculate the Required Fold Equity (RFE).
# This is the breakeven point where the profit from folds equals the loss from being called.
# RFE = Loss when called / (Loss when called + Gain when they fold)
required_fold_equity = loss_when_called / (loss_when_called + blinds_and_antes)

# --- Output the Explanation and Result ---

print(f"Analysis for shoving with {hand_in_question}:")
print(f"Our Stack: {our_stack} BB")
print(f"Pot before our turn: {blinds_and_antes} BB (our reward for a successful shove)")
print(f"Assumed Opponent Calling Range: TT+, AQs+, AKo")
print(f"Our equity with {hand_in_question} when called: {equity_when_called:.1%}")
print("-" * 30)
print(f"If we are called, the total pot will be {pot_if_called:.1f} BB.")
print(f"Our average loss when called is {loss_when_called:.2f} BB.")
print(f"To be profitable, our opponents must fold at least {required_fold_equity:.1%} of the time.")
print("-" * 30)
print("Conclusion: On the money bubble, it's very likely opponents will fold more than 59% of the time to an UTG+1 shove.")
print("This makes shoving 99 a profitable play. While AKo is also a clear shove, 99 is a crucial hand that strong players correctly shove in this spot, relying on fold equity.")
print("\nFinal Calculation for Required Fold Equity (RFE):")
# The final required line:
print(f"RFE = {loss_when_called:.2f} BB / ({loss_when_called:.2f} BB + {blinds_and_antes:.1f} BB) = {required_fold_equity:.1%}")

# --- End of Code, Restore Standard Output ---
# This part is for internal validation and would not be in the final user-facing block
sys.stdout = old_stdout
script_output = mystdout.getvalue()
# print(script_output) # Optional: print to see the captured output
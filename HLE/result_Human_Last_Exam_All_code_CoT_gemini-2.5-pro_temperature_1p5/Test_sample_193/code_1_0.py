import math

# Define the variables from the problem
pot_P = 10
stack_S = 1000 # This is the maximum possible bet size, B

# --- Step 1: Calculate the profit EV of checking with AA ---
# If we check, the most we can expect to win is the current pot.
ev_check_profit = pot_P

# --- Step 2: Calculate the profit EV of betting with AA ---
# To find the optimal bet size, we analyze the EV profit function for betting.
# The function is EV_bet(B) = (2 * P * B) / (P + B).
# This function increases with B, so the optimal bet size is the maximum possible, S.
bet_size_B = stack_S
ev_bet_profit = (2 * pot_P * bet_size_B) / (pot_P + bet_size_B)

# --- Step 3: Determine the optimal strategy ---
# We compare the profit from betting vs. checking.
# If EV of betting is higher, we always bet.
if ev_bet_profit > ev_check_profit:
    # Betting is the superior strategy.
    action = "BET"
    # The sizing is the one that maximizes the profit.
    sizing = bet_size_B
    # Since betting is strictly better, the frequency is 100%.
    frequency = 100
else:
    # In a different scenario, checking or mixing might be optimal.
    # But based on our calculation, this branch will not be taken.
    action = "CHECK"
    sizing = 0
    frequency = 100

# --- Step 4: Format the final answer ---
# Round sizing and frequency to the nearest even number as requested.
sizing_rounded = round(sizing / 2) * 2
frequency_rounded = round(frequency / 2) * 2

# Output the explanation and final answer
print("Calculating the optimal strategy for AA...")
print(f"Pot size (P): {pot_P}")
print(f"Our stack and max bet size (B): {stack_S}")
print("\n--- Comparing Profitability of Actions ---")
print(f"Equation for EV_profit(Check) = P")
print(f"EV_profit when Checking = {ev_check_profit}")
print("\nEquation for EV_profit(Bet) = (2 * P * B) / (P + B)")
print(f"EV_profit when Betting max ({bet_size_B}) = (2 * {pot_P} * {bet_size_B}) / ({pot_P} + {bet_size_B}) = {ev_bet_profit:.2f}")
print(f"\nSince {ev_bet_profit:.2f} > {ev_check_profit}, betting is always the correct action with AA.")
print("The bet size should be maximized to maximize profit.")
print("\nFinal optimal strategy for AA:")
print(f"{action} {sizing_rounded} {frequency_rounded}%")

# The final answer in the required format is on the last line.
# For the purpose of the final output, here it is again.
# The user asked for just one code block and no separate copy-pasting.
final_answer = f"<<<{action} {sizing_rounded} {frequency_rounded}%>>>"

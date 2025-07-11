# 2024 Contribution Limits
# Since you turn 50 in 2024, you are eligible for catch-up contributions for the entire year.
limit_401k_standard = 23000
limit_401k_catch_up = 7500
total_limit_401k = limit_401k_standard + limit_401k_catch_up

limit_ira_standard = 7000
limit_ira_catch_up = 1000
total_limit_ira = limit_ira_standard + limit_ira_catch_up

limit_hsa_individual = 4150

# Your Current Contributions
contribution_401k = 23000
contribution_ira = 5000
contribution_hsa = 4150

# --- Calculations ---

# 1. Calculate remaining 401(k) contribution room
remaining_401k = total_limit_401k - contribution_401k

# 2. Calculate remaining IRA contribution room
remaining_ira = total_limit_ira - contribution_ira

# 3. Calculate remaining HSA contribution room
remaining_hsa = limit_hsa_individual - contribution_hsa

# 4. Calculate total remaining contribution room for retirement accounts
total_remaining = remaining_401k + remaining_ira

# --- Output ---

print("--- 2024 Retirement Contribution Analysis ---")

print(f"\n401(k) Analysis:")
print(f"  Your total 401(k) contribution limit (including ${limit_401k_catch_up} catch-up): ${total_limit_401k}")
print(f"  You have contributed: ${contribution_401k}")
print(f"  Remaining 401(k) contribution room: ${remaining_401k}")

print(f"\nIRA Analysis:")
print(f"  Your total IRA contribution limit (including ${limit_ira_catch_up} catch-up): ${total_limit_ira}")
print(f"  You have contributed: ${contribution_ira}")
print(f"  Remaining IRA contribution room: ${remaining_ira}")

print(f"\nHSA Analysis:")
print(f"  Your HSA is fully funded for 2024 at ${contribution_hsa}.")

print("\n--- Total Additional Contribution Possible ---")
print("The total additional amount you can contribute to your retirement accounts is calculated as:")
print(f"Remaining 401(k) Room + Remaining IRA Room = Total")
print(f"${remaining_401k} + ${remaining_ira} = ${total_remaining}")

print(f"\nYou can legally contribute ${total_remaining} more to your retirement accounts in 2024.")
print("\nFinal Equation:")
print(f"{remaining_401k} + {remaining_ira} = {total_remaining}")
<<<10500>>>
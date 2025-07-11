import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
string_buffer = io.StringIO()
# Redirect stdout to the string buffer
sys.stdout = string_buffer

# User's information
age = 49
birthday_in_2024 = True # Turns 50 in 2024
current_401k_contribution = 23000
current_ira_contribution = 5000
current_hsa_contribution = 4150
current_fsa_contribution = 3200

# 2024 IRS Limits
limit_401k = 23000
catch_up_401k = 7500
limit_ira = 7000
catch_up_ira = 1000
limit_hsa_self_only = 4150
limit_fsa = 3200

# Determine eligibility for catch-up contributions
# The user is eligible if they turn 50 at any point during the calendar year.
eligible_for_catch_up = (age >= 50) or (age == 49 and birthday_in_2024)

# --- 401(k) Calculation ---
total_401k_limit = limit_401k
if eligible_for_catch_up:
    total_401k_limit += catch_up_401k
remaining_401k = max(0, total_401k_limit - current_401k_contribution)

# --- IRA Calculation ---
total_ira_limit = limit_ira
if eligible_for_catch_up:
    total_ira_limit += catch_up_ira
remaining_ira = max(0, total_ira_limit - current_ira_contribution)

# --- Total Calculation ---
total_remaining_contribution = remaining_401k + remaining_ira

# --- Print Explanation ---
print("Based on 2024 IRS rules, here is how much more you can contribute to your retirement accounts this year.")
print("-" * 70)

print("\n1. 401(k) Account:")
print(f"The standard employee contribution limit is ${limit_401k:,}.")
if eligible_for_catch_up:
    print(f"Because you turn 50 in 2024, you are eligible for an additional catch-up contribution of ${catch_up_401k:,}.")
    print(f"Your total personal 401(k) limit is: ${limit_401k:,} (standard) + ${catch_up_401k:,} (catch-up) = ${total_401k_limit:,}.")
else:
    print(f"Your total personal 401(k) limit is ${total_401k_limit:,}.")
print(f"You have contributed ${current_401k_contribution:,} so far.")
print(f"Remaining 401(k) contribution for 2024: ${total_401k_limit:,} - ${current_401k_contribution:,} = ${remaining_401k:,}.")
print("\n(Note: Your employer's match does not count against this personal limit.)")


print("\n2. IRA (Backdoor Roth) Account:")
print(f"The standard IRA contribution limit is ${limit_ira:,}.")
if eligible_for_catch_up:
    print(f"You are also eligible for the IRA catch-up contribution of ${catch_up_ira:,}.")
    print(f"Your total IRA limit is: ${limit_ira:,} (standard) + ${catch_up_ira:,} (catch-up) = ${total_ira_limit:,}.")
else:
    print(f"Your total IRA limit is ${total_ira_limit:,}.")
print(f"You have contributed ${current_ira_contribution:,} so far.")
print(f"Remaining IRA contribution for 2024: ${total_ira_limit:,} - ${current_ira_contribution:,} = ${remaining_ira:,}.")

print("\n3. Other Accounts:")
print(f"You have already contributed the maximum of ${limit_hsa_self_only:,} to your HSA for self-only coverage.")
print(f"You have also contributed the maximum of ${limit_fsa:,} to your FSA.")

print("-" * 70)
print("\nTotal Additional Retirement Contribution for 2024:")
print("This is the sum of the remaining room in your 401(k) and IRA.")
print(f"Total = (Remaining 401(k) room) + (Remaining IRA room)")
print(f"Total = ${remaining_401k:,} + ${remaining_ira:,}")
print(f"You can contribute a total of ${total_remaining_contribution:,} more to your retirement accounts.")

# Get the content from the buffer
output = string_buffer.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Print the content to the actual console
print(output)

# Final answer in specified format
print(f"<<<{total_remaining_contribution}>>>")
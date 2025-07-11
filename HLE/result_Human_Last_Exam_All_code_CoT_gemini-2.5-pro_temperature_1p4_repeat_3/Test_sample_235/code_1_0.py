# 2024 IRS limits for individuals under age 50
ira_limit = 7000
employee_401k_limit = 23000
total_401k_limit = 69000

# User's provided contributions
user_ira_contribution = 5000
user_401k_contribution = 23000
employer_match_rate = 0.50

# --- IRA Calculation ---
# Calculate the remaining amount the user can contribute to an IRA.
additional_ira_allowed = ira_limit - user_ira_contribution

# --- 401(k) Calculation ---
# 1. Calculate the employer match. The problem states a 50% match.
employer_match_amount = user_401k_contribution * employer_match_rate

# 2. Calculate the total current contributions to the 401k (employee + employer).
total_current_401k_contributions = user_401k_contribution + employer_match_amount

# 3. Calculate the remaining space under the overall 401(k) limit ($69,000).
# This space can be filled with after-tax contributions (Mega Backdoor Roth), if the user's plan allows it.
additional_401k_allowed = total_401k_limit - total_current_401k_contributions

# --- Final Calculation ---
# Sum the additional allowed contributions for the final answer.
total_additional_contribution = additional_ira_allowed + additional_401k_allowed

print("This calculation determines the total additional amount you are legally allowed to contribute to your retirement accounts in 2024.")
print("It is the sum of the remaining room in your IRA and the remaining overall space in your 401(k) for after-tax contributions.")
print("\n--- Final Equation ---")
# The final output displays each component of the calculation as requested.
print(f"{int(additional_ira_allowed)} (additional IRA contribution) + {int(additional_401k_allowed)} (additional 401(k) after-tax contribution) = {int(total_additional_contribution)}")
print("\nNote: The ability to make after-tax 401(k) contributions is dependent on your specific employer's plan rules.")
<<<36500>>>
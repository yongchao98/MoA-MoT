import locale

# Set locale to format currency nicely
try:
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
except locale.Error:
    locale.setlocale(locale.LC_ALL, '') # Fallback for other systems

# --- 2024 Contribution Limits (under age 50) ---
employee_401k_limit = 23000
ira_limit = 7000
hsa_limit_self = 4150
fsa_limit = 3200

# --- User's Contributions ---
user_401k_contrib = 23000
user_ira_contrib = 5000
user_hsa_contrib = 4150
user_fsa_contrib = 3200

# --- Calculations ---
remaining_401k = employee_401k_limit - user_401k_contrib
remaining_ira = ira_limit - user_ira_contrib
total_remaining = remaining_401k + remaining_ira

# --- Output ---
print("Here is the analysis of your remaining retirement contributions for 2024:")
print("-" * 50)

# 401k Analysis
print("\n1. 401k Account:")
print(f"The 2024 employee contribution limit for a 401k is {locale.currency(employee_401k_limit, grouping=True)}.")
print(f"You have already contributed {locale.currency(user_401k_contrib, grouping=True)}.")
print(f"Your remaining personal contribution to your 401k is {locale.currency(remaining_401k, grouping=True)}.")
print("(Note: Employer matches do not count against your personal limit.)")
print("-" * 50)


# IRA Analysis
print("\n2. IRA (Traditional/Roth/Backdoor):")
print(f"The 2024 contribution limit for all IRAs is {locale.currency(ira_limit, grouping=True)}.")
print(f"You have contributed {locale.currency(user_ira_contrib, grouping=True)}.")
print("Your remaining personal contribution is calculated as:")
# This line prints each number in the final equation as requested
print(f"{locale.currency(ira_limit, grouping=True)} (Limit) - {locale.currency(user_ira_contrib, grouping=True)} (Contributed) = {locale.currency(remaining_ira, grouping=True)}")
print("-" * 50)


# Other Accounts Analysis
print("\n3. Other Accounts (HSA/FSA):")
print(f"You have contributed {locale.currency(user_hsa_contrib, grouping=True)} to your HSA, which meets the {locale.currency(hsa_limit_self, grouping=True)} maximum for self-only coverage.")
print(f"You have contributed {locale.currency(user_fsa_contrib, grouping=True)} to your FSA, which meets the {locale.currency(fsa_limit, grouping=True)} maximum.")
print("-" * 50)


# Final Summary
print("\n--- Summary ---")
print("Based on this analysis, the total additional amount you are legally allowed to contribute to your retirement accounts for the 2024 tax year is:")
print(f"{locale.currency(total_remaining, grouping=True)}")

import locale

def calculate_retirement_contributions():
    """
    Calculates and explains the remaining retirement contribution room for 2024
    based on the user's provided details.
    """
    # User's information
    user_age_end_of_year = 50
    user_contrib_401k = 23000
    user_contrib_ira = 5000
    user_contrib_hsa = 4150

    # 2024 IRS Contribution Limits
    limit_401k_standard = 23000
    limit_401k_catchup_age_50_plus = 7500
    limit_ira_standard = 7000
    limit_ira_catchup_age_50_plus = 1000
    limit_hsa_self_only = 4150

    # Use locale for formatting numbers with commas
    locale.setlocale(locale.LC_ALL, '')

    # --- 1. 401(k) Contribution Room ---
    # Eligibility for catch-up is determined if you are age 50 by the end of the calendar year.
    total_401k_limit = 0
    if user_age_end_of_year >= 50:
        total_401k_limit = limit_401k_standard + limit_401k_catchup_age_50_plus
    else:
        total_401k_limit = limit_401k_standard
    
    remaining_401k = total_401k_limit - user_contrib_401k

    # --- 2. IRA Contribution Room ---
    total_ira_limit = 0
    if user_age_end_of_year >= 50:
        total_ira_limit = limit_ira_standard + limit_ira_catchup_age_50_plus
    else:
        total_ira_limit = limit_ira_standard

    remaining_ira = total_ira_limit - user_contrib_ira

    # --- 3. HSA & FSA ---
    remaining_hsa = limit_hsa_self_only - user_contrib_hsa
    remaining_hsa = max(0, remaining_hsa) # Cannot be negative

    # --- 4. Total Remaining Contribution ---
    total_remaining = remaining_401k + remaining_ira

    # --- Output ---
    print("Based on the 2024 tax year rules, here is a breakdown of your remaining retirement contribution room:")
    
    print("\n--- 401(k) ---")
    print(f"The standard employee contribution limit is ${locale.format_string('%d', limit_401k_standard, grouping=True)}.")
    print(f"Because you turn 50 in 2024, you are eligible for an additional catch-up contribution of ${locale.format_string('%d', limit_401k_catchup_age_50_plus, grouping=True)}.")
    print(f"Your total 401(k) contribution limit is ${locale.format_string('%d', total_401k_limit, grouping=True)}.")
    print(f"You have already contributed ${locale.format_string('%d', user_contrib_401k, grouping=True)}.")
    print(f"Remaining 401(k) contribution allowed: ${locale.format_string('%d', remaining_401k, grouping=True)}")

    print("\n--- IRA (Backdoor Roth) ---")
    print(f"The standard IRA contribution limit is ${locale.format_string('%d', limit_ira_standard, grouping=True)}.")
    print(f"You are also eligible for the age 50+ catch-up contribution of ${locale.format_string('%d', limit_ira_catchup_age_50_plus, grouping=True)}.")
    print(f"Your total IRA contribution limit is ${locale.format_string('%d', total_ira_limit, grouping=True)}.")
    print(f"You have already contributed ${locale.format_string('%d', user_contrib_ira, grouping=True)}.")
    print(f"Remaining IRA contribution allowed: ${locale.format_string('%d', remaining_ira, grouping=True)}")

    print("\n--- HSA & FSA ---")
    print(f"You have contributed ${locale.format_string('%d', user_contrib_hsa, grouping=True)} to your HSA, which is the maximum for self-only coverage (${locale.format_string('%d', limit_hsa_self_only, grouping=True)}) for those under 55. You have ${locale.format_string('%d', remaining_hsa, grouping=True)} remaining to contribute here.")
    print("Your FSA contributions do not count towards these retirement savings limits.")

    print("\n--- Total Additional Contribution for 2024 ---")
    print("To calculate the total additional amount you can contribute to your retirement accounts, we sum the remaining room in your 401(k) and IRA.")
    print(f"Final Calculation: {remaining_401k} (from 401k) + {remaining_ira} (from IRA) = {total_remaining}")

if __name__ == '__main__':
    calculate_retirement_contributions()
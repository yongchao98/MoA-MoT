def calculate_remaining_retirement_contributions():
    """
    Calculates and displays the remaining retirement contribution room for a user in 2024.
    """
    # --- 2024 IRS Contribution Limits ---
    LIMIT_401K_REGULAR = 23000
    LIMIT_401K_CATCH_UP_50_PLUS = 7500
    LIMIT_IRA_REGULAR = 7000
    LIMIT_IRA_CATCH_UP_50_PLUS = 1000
    LIMIT_HSA_SELF_ONLY = 4150
    LIMIT_FSA = 3200

    # --- User's Provided Information ---
    user_contribution_401k = 23000
    user_contribution_ira = 5000
    user_contribution_hsa = 4150
    user_contribution_fsa = 3200

    # Since you turn 50 in 2024, you are eligible for catch-up contributions.
    max_401k_limit = LIMIT_401K_REGULAR + LIMIT_401K_CATCH_UP_50_PLUS
    max_ira_limit = LIMIT_IRA_REGULAR + LIMIT_IRA_CATCH_UP_50_PLUS

    # --- Calculate Remaining Contribution Room ---
    remaining_401k = max_401k_limit - user_contribution_401k
    remaining_ira = max_ira_limit - user_contribution_ira
    total_remaining_retirement = remaining_401k + remaining_ira

    # --- Print the Detailed Breakdown ---
    print("Here is a breakdown of your remaining 2024 retirement contribution allowances.")
    print("-" * 70)
    print("First, let's review your maximums based on your age:")
    print("Since you turn 50 in 2024, you are eligible for 'catch-up' contributions.")
    
    print(f"\nYour Maximum 401(k) Employee Contribution: ${LIMIT_401K_REGULAR} (standard) + ${LIMIT_401K_CATCH_UP_50_PLUS} (catch-up) = ${max_401k_limit}")
    print(f"Your Maximum IRA Contribution: ${LIMIT_IRA_REGULAR} (standard) + ${LIMIT_IRA_CATCH_UP_50_PLUS} (catch-up) = ${max_ira_limit}")
    print("-" * 70)
    
    print("Now, let's calculate how much more you can legally contribute:\n")

    # 401(k) Calculation
    print("401(k) Account:")
    print(f"  Remaining contribution allowed = ${max_401k_limit} (Your Limit) - ${user_contribution_401k} (Contributed) = ${remaining_401k}")

    # IRA Calculation
    print("\nIRA Account (via Backdoor Roth):")
    print(f"  Remaining contribution allowed = ${max_ira_limit} (Your Limit) - ${user_contribution_ira} (Contributed) = ${remaining_ira}")

    # Notes on other accounts
    print("\nOther Tax-Advantaged Accounts:")
    print(f"  HSA & FSA: You have already contributed the maximum amounts of ${user_contribution_hsa} and ${user_contribution_fsa} respectively for 2024.")
    print("-" * 70)

    # Final Total Calculation
    print("Total additional amount you can contribute to your retirement accounts in 2024:")
    print(f"${remaining_401k} (from 401k) + ${remaining_ira} (from IRA) = ${total_remaining_retirement}")

# Run the calculation
calculate_remaining_retirement_contributions()
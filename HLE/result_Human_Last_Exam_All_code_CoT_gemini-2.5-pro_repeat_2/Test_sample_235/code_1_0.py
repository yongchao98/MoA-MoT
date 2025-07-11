def calculate_retirement_contributions():
    """
    Calculates the remaining retirement contribution room for a 49-year-old for the 2024 tax year.
    """
    # User's Information
    # Although the user is 49, they turn 50 in 2024, making them eligible for catch-up contributions.
    user_is_eligible_for_50_plus_catchup = True
    user_contribution_401k = 23000
    user_contribution_ira = 5000
    user_contribution_hsa = 4150

    # 2024 IRS Limits
    limit_401k_employee = 23000
    catchup_401k = 7500
    limit_ira = 7000
    catchup_ira = 1000
    limit_hsa_individual = 4150

    # --- 401(k) Calculation ---
    total_limit_401k = limit_401k_employee
    if user_is_eligible_for_50_plus_catchup:
        total_limit_401k += catchup_401k
    remaining_401k = max(0, total_limit_401k - user_contribution_401k)

    # --- IRA Calculation ---
    total_limit_ira = limit_ira
    if user_is_eligible_for_50_plus_catchup:
        total_limit_ira += catchup_ira
    remaining_ira = max(0, total_limit_ira - user_contribution_ira)

    # --- HSA Calculation ---
    # HSA catch-up is for age 55+, so it does not apply here.
    total_limit_hsa = limit_hsa_individual
    remaining_hsa = max(0, total_limit_hsa - user_contribution_hsa)

    # --- Total Calculation ---
    total_remaining_contribution = remaining_401k + remaining_ira + remaining_hsa

    # --- Print the explanation ---
    print("Here is your remaining retirement contribution breakdown for 2024:\n")

    print("--- 401(k) ---")
    print(f"Your total 401(k) limit (including a ${catchup_401k} catch-up): ${total_limit_401k}")
    print(f"You have contributed: ${user_contribution_401k}")
    print(f"Remaining 401(k) contribution room: ${remaining_401k}\n")

    print("--- IRA (Roth/Traditional) ---")
    print(f"Your total IRA limit (including a ${catchup_ira} catch-up): ${total_limit_ira}")
    print(f"You have contributed: ${user_contribution_ira}")
    print(f"Remaining IRA contribution room: ${remaining_ira}\n")

    print("--- HSA (Health Savings Account) ---")
    print(f"Your total HSA limit: ${total_limit_hsa}")
    print(f"You have contributed: ${user_contribution_hsa}")
    print(f"Remaining HSA contribution room: ${remaining_hsa}\n")
    
    print("-" * 40)
    print("Total Remaining Contribution Calculation:")
    print(f"The total amount you can still contribute to these tax-advantaged accounts is:")
    print(f"${remaining_401k} (401k) + ${remaining_ira} (IRA) + ${remaining_hsa} (HSA) = ${total_remaining_contribution}")
    print("-" * 40)

if __name__ == '__main__':
    calculate_retirement_contributions()
def calculate_remaining_retirement_contributions():
    """
    Calculates and explains the remaining retirement contributions for a 49-year-old in 2024.
    """
    # User's information
    age = 49
    # The user turns 50 in 2024, making them eligible for catch-up contributions for the entire year.
    eligible_for_catch_up = True
    
    current_401k_contribution = 23000
    current_ira_contribution = 5000
    current_hsa_contribution = 4150
    current_fsa_contribution = 3200

    # 2024 IRS Contribution Limits
    LIMIT_401K_STANDARD = 23000
    LIMIT_401K_CATCH_UP = 7500
    LIMIT_IRA_STANDARD = 7000
    LIMIT_IRA_CATCH_UP = 1000
    LIMIT_HSA_SELF = 4150
    LIMIT_FSA = 3200

    # --- Calculations ---

    # 1. Calculate total 401(k) limit for the user
    total_401k_limit = LIMIT_401K_STANDARD
    if eligible_for_catch_up:
        total_401k_limit += LIMIT_401K_CATCH_UP
    
    # 2. Calculate remaining 401(k) contribution room
    remaining_401k = total_401k_limit - current_401k_contribution

    # 3. Calculate total IRA limit for the user
    total_ira_limit = LIMIT_IRA_STANDARD
    if eligible_for_catch_up:
        total_ira_limit += LIMIT_IRA_CATCH_UP

    # 4. Calculate remaining IRA contribution room
    remaining_ira = total_ira_limit - current_ira_contribution
    
    # 5. Calculate total remaining retirement contributions
    total_remaining = remaining_401k + remaining_ira

    # --- Output ---
    print("--- 2024 Retirement Contribution Analysis ---")
    
    print(f"\nYour Total 401(k) Contribution Limit (including ${LIMIT_401K_CATCH_UP} catch-up): ${total_401k_limit}")
    print(f"You have contributed: ${current_401k_contribution}")
    print(f"Remaining 401(k) contribution room: ${remaining_401k}")

    print(f"\nYour Total IRA Contribution Limit (including ${LIMIT_IRA_CATCH_UP} catch-up): ${total_ira_limit}")
    print(f"You have contributed: ${current_ira_contribution}")
    print(f"Remaining IRA contribution room: ${remaining_ira}")
    
    print(f"\nNote: Your HSA is fully funded at ${current_hsa_contribution}/${LIMIT_HSA_SELF}, and your FSA is not a retirement account.")
    
    print("\n--- Final Calculation ---")
    print("Total additional money you can legally contribute to your retirement accounts in 2024:")
    # Final equation format as requested
    print(f"${remaining_401k} (401k) + ${remaining_ira} (IRA) = ${total_remaining}")

calculate_remaining_retirement_contributions()
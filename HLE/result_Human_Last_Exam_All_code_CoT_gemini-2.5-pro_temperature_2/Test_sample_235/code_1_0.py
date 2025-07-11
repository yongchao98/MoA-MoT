def calculate_remaining_retirement_contributions():
    """
    Calculates the remaining retirement contribution room for a user in 2024.
    """
    # 2024 Contribution Limits
    LIMIT_401K_STANDARD = 23000
    LIMIT_401K_CATCHUP_50_PLUS = 7500
    LIMIT_IRA_STANDARD = 7000
    LIMIT_IRA_CATCHUP_50_PLUS = 1000

    # User's Contributions
    user_contrib_401k = 23000
    user_contrib_ira = 5000

    # The user turns 50 in 2024, so they are eligible for catch-up contributions for the entire year.
    total_401k_limit = LIMIT_401K_STANDARD + LIMIT_401K_CATCHUP_50_PLUS
    total_ira_limit = LIMIT_IRA_STANDARD + LIMIT_IRA_CATCHUP_50_PLUS

    # Calculate remaining contribution room for each account
    remaining_401k = total_401k_limit - user_contrib_401k
    remaining_ira = total_ira_limit - user_contrib_ira
    total_remaining = remaining_401k + remaining_ira

    # Print the explanation
    print("Based on the 2024 tax year limits, here is your remaining contribution room:\n")
    print("--- 401(k) Calculation ---")
    print(f"Your total 401(k) contribution limit (age 50+) is ${total_401k_limit:,} (${LIMIT_401K_STANDARD:,} standard + ${LIMIT_401K_CATCHUP_50_PLUS:,} catch-up).")
    print(f"You have already contributed ${user_contrib_401k:,}.")
    print(f"Remaining 401(k) Contribution Room: ${remaining_401k:,}\n")

    print("--- IRA Calculation ---")
    print(f"Your total IRA contribution limit (age 50+) is ${total_ira_limit:,} (${LIMIT_IRA_STANDARD:,} standard + ${LIMIT_IRA_CATCHUP_50_PLUS:,} catch-up).")
    print(f"You have already contributed ${user_contrib_ira:,} to your backdoor Roth IRA.")
    print(f"Remaining IRA Contribution Room: ${remaining_ira:,}\n")

    print("--- Total Remaining Contribution ---")
    print(f"You can legally contribute an additional ${remaining_401k:,} to your 401(k) and ${remaining_ira:,} to your IRA.")
    print(f"Total additional contribution allowed = ${remaining_401k:,} (401k) + ${remaining_ira:,} (IRA) = ${total_remaining:,}")

    # Note on other accounts
    print("\nNote: Your HSA is fully funded for 2024, and your FSA is not classified as a retirement account.")
    
    # Return final answer for the system
    # Do not manually print this final line, the function call will do it.
    return total_remaining

# Execute the function and print the final formatted answer
final_answer = calculate_remaining_retirement_contributions()
print(f"\n<<<{final_answer}>>>")

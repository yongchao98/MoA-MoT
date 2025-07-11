def calculate_remaining_retirement_contributions():
    """
    Calculates and explains the remaining retirement contribution room for the user in 2024.
    """
    # 2024 IRS Limits
    LIMIT_401K_EMPLOYEE = 23000
    LIMIT_401K_OVERALL = 69000
    LIMIT_IRA = 7000
    LIMIT_HSA_SELF = 4150

    # User's current contributions
    user_contribution_401k = 23000
    user_contribution_ira = 5000
    user_contribution_hsa = 4150
    user_contribution_fsa = 3200

    # User's employer match details
    employer_match_rate = 0.50

    # --- Calculations ---

    # 1. IRA Calculation
    remaining_ira_contribution = LIMIT_IRA - user_contribution_ira

    # 2. 401k Calculation
    # User has maxed out the standard employee contribution
    remaining_employee_401k_contribution = LIMIT_401K_EMPLOYEE - user_contribution_401k

    # Calculate employer match
    employer_match_contribution = user_contribution_401k * employer_match_rate

    # Calculate total current 401k contributions against the overall limit
    total_current_401k_contributions = user_contribution_401k + employer_match_contribution

    # Calculate remaining room in the overall 401k limit (for after-tax contributions)
    remaining_overall_401k_contribution = LIMIT_401K_OVERALL - total_current_401k_contributions
    
    # 3. Total Calculation
    total_remaining_contribution = remaining_ira_contribution + remaining_overall_401k_contribution

    # --- Print Explanation ---
    print("Based on the 2024 IRS limits, here is your remaining contribution allowance:")
    print("-" * 60)
    
    print(f"\n1. Individual Retirement Account (IRA) Analysis:")
    print(f"The annual IRA contribution limit is ${LIMIT_IRA:,}.")
    print(f"You have contributed ${user_contribution_ira:,}.")
    print(f"Remaining IRA Contribution = ${LIMIT_IRA:,} (Limit) - ${user_contribution_ira:,} (Your Contribution)")
    print(f"You can contribute an additional ${remaining_ira_contribution:,.2f} to your IRA.")

    print("\n2. 401k Plan Analysis:")
    print("You have already contributed the maximum employee amount of $23,000 for 2024.")
    print("However, we can calculate additional space under the overall plan limit of $69,000.")
    
    print("\n   Calculating Your Total 401k Contributions So Far:")
    print(f"   Your Employee Contribution: ${user_contribution_401k:,.2f}")
    print(f"   Your Employer Match (50% of your contribution): ${employer_match_contribution:,.2f}")
    print(f"   Total Current Contributions = ${user_contribution_401k:,.2f} + ${employer_match_contribution:,.2f} = ${total_current_401k_contributions:,.2f}")

    print("\n   Calculating Remaining After-Tax 401k Space:")
    print(f"   Overall 401k Limit (Employee + Employer + After-Tax): ${LIMIT_401K_OVERALL:,.2f}")
    print(f"   Remaining Space = ${LIMIT_401K_OVERALL:,.2f} (Limit) - ${total_current_401k_contributions:,.2f} (Total Current)")
    print(f"   You can contribute an additional ${remaining_overall_401k_contribution:,.2f} to your 401k as after-tax contributions.")

    print("-" * 60)
    print("\n3. Final Calculation:")
    print("Your total additional retirement contribution capacity for 2024 is the sum of the remaining IRA and 401k allowances.")
    print(f"\nTotal = ${remaining_ira_contribution:,.2f} (Remaining IRA) + ${remaining_overall_401k_contribution:,.2f} (Remaining 401k)")
    print(f"\nYou can legally contribute an additional ${total_remaining_contribution:,.2f} to your retirement accounts in 2024.")

    print("\nNote: Your HSA is fully funded for the year, and your FSA is not a retirement account. This calculation assumes your 401k plan allows for after-tax contributions.")
    
    # Return value for the final answer format
    return total_remaining_contribution

# Run the calculation and store the result
final_answer = calculate_remaining_retirement_contributions()
print(f"\n<<<{final_answer}>>>")

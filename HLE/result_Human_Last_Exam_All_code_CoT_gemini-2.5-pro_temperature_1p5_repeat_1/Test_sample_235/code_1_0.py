def calculate_remaining_retirement_contributions():
    """
    Calculates the remaining legally allowed retirement contributions for the user in 2024.
    """
    # --- 2024 IRS Contribution Limits (for individuals under 50) ---
    ira_limit = 7000
    k401_employee_limit = 23000
    k401_overall_limit = 69000  # Employee + Employer + Other
    hsa_limit_single = 4150

    # --- User's Provided Contributions ---
    user_k401_contribution = 23000
    user_ira_contribution = 5000
    user_hsa_contribution = 4150
    user_fsa_contribution = 3200 # Note: FSA is not a retirement account

    print("Analyzing your 2024 retirement contribution limits...\n")

    # --- 1. IRA Contribution Analysis ---
    print("--- IRA Contributions ---")
    remaining_ira_room = ira_limit - user_ira_contribution
    print(f"The 2024 IRA contribution limit for someone under 50 is ${ira_limit:,.2f}.")
    print(f"You have contributed ${user_ira_contribution:,.2f}.")
    if remaining_ira_room > 0:
        print(f"You can contribute an additional ${remaining_ira_room:,.2f} to your IRA.")
    else:
        print("You have maxed out your IRA contributions for the year.")
    print("-" * 25 + "\n")


    # --- 2. HSA Contribution Analysis ---
    print("--- HSA Contributions ---")
    remaining_hsa_room = hsa_limit_single - user_hsa_contribution
    print(f"The 2024 HSA contribution limit for self-only coverage is ${hsa_limit_single:,.2f}.")
    print(f"You have contributed ${user_hsa_contribution:,.2f}.")
    if remaining_hsa_room <= 0:
      print("You have maxed out your HSA contributions for the year.")
    print("-" * 25 + "\n")


    # --- 3. 401(k) Contribution Analysis ---
    print("--- 401(k) Contributions ---")
    # Employee Deferral Limit
    remaining_employee_k401_room = k401_employee_limit - user_k401_contribution
    print(f"The 2024 employee 401(k) contribution limit is ${k401_employee_limit:,.2f}.")
    print(f"You have contributed ${user_k401_contribution:,.2f}, so you have maxed out your standard employee contributions.")
    
    print("\nHowever, we can check for additional room under the overall 401(k) limit.")
    print(f"The total 401(k) limit from all sources (employee, employer, after-tax) is ${k401_overall_limit:,.2f}.")
    
    # Employer Match Calculation
    # Assumption: The employer matches 50% of the employee's contribution.
    employer_match_rate = 0.50
    employer_match = user_k401_contribution * employer_match_rate
    print("\nCalculating your employer match...")
    print(f"Based on your contribution of ${user_k401_contribution:,.2f} and a 50% match rate, your employer contributed ${employer_match:,.2f}.")
    print("(Note: This assumes your employer matches 50% of your total contribution. Check your plan documents to confirm the exact match formula.)")

    # Overall 401(k) Room Calculation (for "Mega Backdoor Roth")
    total_current_k401_contributions = user_k401_contribution + employer_match
    remaining_after_tax_k401_room = k401_overall_limit - total_current_k401_contributions
    
    print(f"\nYour total 401(k) contributions so far are: ${user_k401_contribution:,.2f} (you) + ${employer_match:,.2f} (employer) = ${total_current_k401_contributions:,.2f}.")
    if remaining_after_tax_k401_room > 0:
        print(f"This leaves ${remaining_after_tax_k401_room:,.2f} of potential contribution room under the overall limit.")
        print("If your plan allows, this amount can be contributed as after-tax dollars (often called a 'Mega Backdoor Roth' strategy).")
    else:
        print("Your total contributions have reached the overall limit for the year.")
    print("-" * 25 + "\n")

    # --- 4. Final Summary ---
    print("--- Total Remaining Contribution Summary ---")
    total_remaining_contribution = remaining_ira_room + remaining_after_tax_k401_room
    print("To calculate your total remaining contribution space, we sum the available room from each account:")
    print(f"Remaining IRA Room + Potential After-Tax 401(k) Room = Total")
    print(f"${remaining_ira_room:,.2f} + ${remaining_after_tax_k401_room:,.2f} = ${total_remaining_contribution:,.2f}")
    
    print("\n*Final Answer:*")
    print(f"You can legally contribute a total of ${total_remaining_contribution:,.2f} more to your retirement accounts in 2024.")
    print("This consists of $2,000.00 to your IRA and, if your 401(k) plan allows it, $34,500.00 in after-tax contributions to your 401(k).")

    return total_remaining_contribution

# Execute the function to get the final result.
final_answer = calculate_remaining_retirement_contributions()
# The final answer is wrapped for the system. The printed output above provides the detailed explanation.
# print(f"<<<{final_answer}>>>")
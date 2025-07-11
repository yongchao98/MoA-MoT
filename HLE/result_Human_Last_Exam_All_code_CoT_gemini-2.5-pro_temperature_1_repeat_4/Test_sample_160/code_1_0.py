import datetime

def solve_bank_guarantee_case():
    """
    Analyzes a bank guarantee case based on URDG 758 rules
    to determine the correct course of action for the presenting bank.
    """
    # --- Data from the problem ---
    presentation_date_str = "2020-04-27"
    holidays_str = ["2020-04-30", "2020-05-01"]
    refusal_date_str = "2020-05-06"
    
    # --- Convert strings to date objects ---
    presentation_date = datetime.datetime.strptime(presentation_date_str, "%Y-%m-%d").date()
    holidays = [datetime.datetime.strptime(d, "%Y-%m-%d").date() for d in holidays_str]
    refusal_date = datetime.datetime.strptime(refusal_date_str, "%Y-%m-%d").date()

    # --- Step 1: Calculate the examination deadline ---
    print("Step 1: Calculate the 5-business-day examination deadline for the issuing bank.")
    print("-----------------------------------------------------------------------------")
    print(f"Presentation of documents was on: {presentation_date.strftime('%A, %d %B %Y')}")
    print("According to URDG 758, the examination period of 5 business days starts on the day following the day of presentation.")
    
    business_days_count = 0
    current_date = presentation_date
    examination_days = []

    # Loop to find the 5 business days for examination
    while business_days_count < 5:
        current_date += datetime.timedelta(days=1)
        # Weekday check: Monday is 0, Sunday is 6. Skip Saturday (5) and Sunday (6).
        if current_date.weekday() >= 5:
            continue
        # Holiday check
        if current_date in holidays:
            continue
        
        business_days_count += 1
        examination_days.append(current_date)
        
    deadline_date = examination_days[-1]
    
    print("\nCalculation of the 5 business days for examination:")
    for i, day in enumerate(examination_days, 1):
        print(f"Business Day {i}: {day.strftime('%A, %d %B %Y')}")

    print(f"\nConclusion 1: The deadline for the issuing bank to send a refusal was the close of business on {deadline_date.strftime('%A, %d %B %Y')}.")
    print(f"The bank sent its refusal messages on {refusal_date.strftime('%A, %d %B %Y')}, which is within the 5-business-day limit.\n")

    # --- Step 2: Analyze the validity of the refusal messages ---
    print("Step 2: Analyze the validity of sending two separate refusal messages.")
    print("--------------------------------------------------------------------")
    print("The governing rule is URDG 758, Article 24 (Notice of Rejection).")
    print("Key points from Article 24:")
    print("  - A notice of rejection must state *each discrepancy* in respect of which the guarantor rejects the demand.")
    print("  - If a notice of rejection fails to state a particular discrepancy, the guarantor is precluded (prevented) from claiming that discrepancy later.")
    
    print("\nIssuing Bank's actions on 06 May 2020:")
    print("1. At 09:00 am: Sent refusal notice for 'wrong demand amount'.")
    print("2. At 03:00 pm: Sent a second refusal notice for 'missing beneficiary’s signature'.")

    print("\nAnalysis:")
    print("By sending the first notice at 09:00 am, the bank officially stated its reasons for refusal. Under URDG 758, this notice must contain all discrepancies.")
    print("The bank's attempt to add a new discrepancy ('missing signature') in a second message is invalid. By issuing the first notice, it waived its right to claim any other discrepancies.")
    
    print("\nConclusion 2: The only legally valid discrepancy that the presenting bank must address is the first one: 'wrong demand amount'.\n")

    # --- Step 3: Determine the best course of action ---
    print("Step 3: Determine the best course of action for the presenting bank.")
    print("------------------------------------------------------------------")
    print("The primary issue to resolve is the 'wrong demand amount'. The second discrepancy should be challenged as invalidly raised.")
    print("Due to the bank's closure (force majeure), the guarantee's validity was extended for 30 days from the reopening date (27 April 2020), so there is still time to submit a corrected presentation.")
    print("The most logical action is to focus on fixing the validly raised issue to secure payment for the beneficiary.")
    
    print("\nRecommended Action: The presenting bank should inform the beneficiary about the first, valid discrepancy ('wrong demand amount') so the beneficiary can make the necessary correction. The corrected demand can then be re-presented.")
    print("\nThis action corresponds to Answer Choice B.")

solve_bank_guarantee_case()
<<<B>>>
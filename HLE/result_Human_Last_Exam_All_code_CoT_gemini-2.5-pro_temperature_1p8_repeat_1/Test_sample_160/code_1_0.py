import datetime

def solve_bank_guarantee_case():
    """
    Analyzes a bank guarantee case based on URDG 758 rules
    to determine the correct course of action for the presenting bank.
    """

    # --- Step 1 & 2: Define dates and calculate the examination period deadline ---
    presentation_date = datetime.date(2020, 4, 27)
    holidays = [datetime.date(2020, 4, 30), datetime.date(2020, 5, 1)] # Vietnamese holiday and Labor Day

    print("--- Case Analysis based on URDG 758 ---")
    print(f"Presentation Date: {presentation_date.strftime('%d %B %Y')} (Monday)")

    business_days_to_examine = 5
    business_days_counted = 0
    current_date = presentation_date

    print("\n--- Calculating the 5-Business-Day Examination Period ---")
    # The examination period starts on the business day following the day of presentation.
    while business_days_counted < business_days_to_examine:
        current_date += datetime.timedelta(days=1)
        # weekday() returns 5 for Saturday and 6 for Sunday
        if current_date.weekday() < 5 and current_date not in holidays:
            business_days_counted += 1
            print(f"Business Day {business_days_counted}: {current_date.strftime('%d %B %Y')}")
        else:
            print(f"Non-Business Day: {current_date.strftime('%d %B %Y')} ({current_date.strftime('%A')})")

    examination_deadline = current_date
    print(f"\n[Equation] Deadline for refusal = Presentation Date + 5 Business Days")
    print(f"Final Deadline for Issuing Bank to Refuse: Close of business on {examination_deadline.strftime('%d %B %Y')}")

    # --- Step 3: Analyze the refusal notices ---
    refusal_date = datetime.date(2020, 5, 6)
    print("\n--- Analyzing the Refusal Notices ---")
    print(f"Refusal messages sent on: {refusal_date.strftime('%d %B %Y')}")

    if refusal_date <= examination_deadline:
        print("Finding 1: The refusal notice was sent within the 5-business-day time limit.")
    else:
        print("Finding 1: The refusal notice was sent after the deadline. The issuing bank must honor the demand.")

    print("\nFinding 2: URDG 758 Article 24(d) requires a single notice of refusal to state ALL discrepancies.")
    print("  - The issuing bank sent a first notice listing 'wrong demand amount'.")
    print("  - By doing so, it waived its right to claim any other discrepancies.")
    print("  - The second notice listing 'missing beneficiary's signature' is therefore invalid.")

    # --- Step 4: Determine the correct action ---
    print("\n--- Conclusion ---")
    print("The only valid discrepancy is 'wrong demand amount' from the first notice.")
    print("The presenting bank's correct action is to inform the beneficiary of this valid discrepancy so it can be corrected.")
    print("\nThis corresponds to Answer Choice B.")

solve_bank_guarantee_case()
<<<B>>>
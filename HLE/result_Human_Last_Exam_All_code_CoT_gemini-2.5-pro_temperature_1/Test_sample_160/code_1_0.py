import datetime

def analyze_trade_finance_scenario():
    """
    Analyzes the demand guarantee scenario based on URDG 758 rules.
    This function calculates the examination period and evaluates the validity
    of the issuing bank's refusal notices to determine the correct course of action.
    """

    # --- Step 1: Define Key Dates and Parameters ---
    presentation_date = datetime.date(2020, 4, 27)
    refusal_date = datetime.date(2020, 5, 6)
    
    # Non-banking days in Vietnam
    # Note: May 1st is also a holiday (International Labour Day)
    holidays = {
        datetime.date(2020, 4, 30), # Vietnamese Reunification Day holiday
        datetime.date(2020, 5, 1)   # International Labour Day
    }
    # Monday is 0, Sunday is 6
    weekend_days = {5, 6} # Saturday, Sunday

    # URDG 758 Rules
    examination_period_days = 5
    
    # --- Step 2: Calculate the Examination Period Deadline ---
    print("--- Analysis of the Bank Guarantee Scenario ---")
    print(f"\n1. Calculating the Examination Period (as per URDG 758 Article 20)")
    print(f"Presentation of documents was on: Monday, {presentation_date.strftime('%d %b %Y')}")
    print(f"The issuing bank has {examination_period_days} business days to examine the documents.")

    business_days_passed = 0
    current_date = presentation_date
    print("Examination days count:")
    while business_days_passed < examination_period_days:
        current_date += datetime.timedelta(days=1)
        day_name = current_date.strftime('%A')
        if current_date.weekday() in weekend_days:
            print(f"- {current_date.strftime('%d %b %Y')} ({day_name}) - Weekend (not counted)")
        elif current_date in holidays:
            print(f"- {current_date.strftime('%d %b %Y')} ({day_name}) - Holiday (not counted)")
        else:
            business_days_passed += 1
            print(f"- {current_date.strftime('%d %b %Y')} ({day_name}) - Business Day #{business_days_passed}")
    
    deadline_date = current_date
    print(f"\nThe deadline for refusal was the close of business on: {deadline_date.strftime('%A, %d %b %Y')}.")

    # --- Step 3: Analyze the Refusal Messages ---
    print("\n2. Analyzing the Refusal (as per URDG 758 Article 24)")
    print(f"The issuing bank sent refusal messages on {refusal_date.strftime('%d %b %Y')}, which is within the deadline.")
    print("However, URDG 758 requires a *single notice* of rejection stating *all* discrepancies.")
    print("The issuing bank sent two separate notices:")
    print("  - 09:00 am: Notice 1 (Wrong demand amount)")
    print("  - 03:00 pm: Notice 2 (Missing beneficiary's signature)")
    print("This 'piecemeal' rejection is a procedural failure and violates the rules.")

    # --- Step 4: Determine the Consequence ---
    print("\n3. Determining the Consequence (as per URDG 758 Article 24(d))")
    print("The consequence of a defective refusal is 'preclusion'.")
    print("This means the issuing bank is 'precluded from claiming that the demand... [is] not in conformity'.")
    print("In short, the issuing bank has lost its right to reject the documents and is now obligated to pay.")

    # --- Step 5: Select the Best Action ---
    print("\n4. Conclusion and Recommended Action")
    print("The presenting bank is in a strong position due to the issuing bank's error.")
    print("The correct action is to assert that the issuing bank is precluded from refusal and must honor the demand.")
    print("Therefore, the presenting bank should ask the issuing bank to honor the documents.")

# Execute the analysis function
analyze_trade_finance_scenario()
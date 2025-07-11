import datetime

def solve_banking_case():
    """
    Analyzes a bank guarantee case to determine the correct course of action.
    """
    # --- Part 1: Define Key Dates and Parameters ---
    presentation_date = datetime.date(2020, 4, 27)
    rejection_1_date = datetime.date(2020, 5, 6)
    examination_period_days = 5
    
    holidays = [
        datetime.date(2020, 4, 30),  # Vietnamese Reunification Day
        datetime.date(2020, 5, 1)   # International Labor Day
    ]

    print("--- Case Analysis ---")
    
    # --- Part 2: Calculate the Refusal Deadline (based on URDG 758 Article 20) ---
    print(f"\nStep 1: Calculating the deadline for the issuing bank to refuse the demand.")
    print(f"Presentation was on Monday, 27 April 2020.")
    print(f"The bank has {examination_period_days} business days following presentation to examine the documents.")
    
    business_days_counted = 0
    current_date = presentation_date
    
    while business_days_counted < examination_period_days:
        current_date += datetime.timedelta(days=1)
        # weekday() is 0 for Monday, 5 for Saturday, 6 for Sunday
        is_weekend = current_date.weekday() >= 5
        is_holiday = current_date in holidays
        
        if not is_weekend and not is_holiday:
            business_days_counted += 1

    rejection_deadline = current_date
    
    print(f"The deadline for refusal was the close of business on: {rejection_deadline.strftime('%A, %d %B %Y')}.")
    print(f"The bank sent its refusal on {rejection_1_date.strftime('%A, %d %B %Y')}, which is within this period.\n")
    
    # --- Part 3: Analyze the Validity of the Refusal (based on URDG 758 Article 24) ---
    print("Step 2: Analyzing the procedure of the refusal.")
    print("According to URDG 758 Article 24, a refusal notice must state ALL discrepancies in a SINGLE notice.")
    
    print("\nFacts of the case:")
    print("  - Notice 1 (09:00 am): Stated discrepancy of 'wrong demand amount'.")
    print("  - Notice 2 (03:00 pm): Stated discrepancy of 'missing beneficiary’s signature'.")
    
    print("\nAnalysis:")
    print("The issuing bank sent two separate notices. This is a clear violation of the 'single notice' rule.")
    
    # --- Part 4: Determine the Consequence and Final Answer ---
    print("\nStep 3: Determining the consequence of the procedural violation.")
    print("According to URDG 758 Article 24(f), if a bank fails to follow the refusal rules, it is 'precluded from claiming that the demand...[is] not in conformity'.")
    
    print("\n--- Final Conclusion ---")
    print("The issuing bank's refusal is invalid due to the procedural failure (sending two notices).")
    print("Therefore, the issuing bank is obligated to pay the demand.")
    print("The presenting bank's correct course of action is to ask the issuing bank to honor the documents.")

solve_banking_case()
<<<A>>>
import datetime

def solve_banking_case():
    """
    Analyzes the bank guarantee case according to URDG 758 rules.
    """
    # --- Step 1: Define key dates and information ---
    presentation_date = datetime.date(2020, 4, 27)
    reopening_date = datetime.date(2020, 4, 27)
    first_refusal_date = datetime.date(2020, 5, 6)
    
    # Non-banking days in Vietnam for the period
    # Saturday and Sunday are non-banking days
    # April 30th is a public holiday (Reunification Day)
    # May 1st is a public holiday (Labour Day)
    holidays = [datetime.date(2020, 4, 30), datetime.date(2020, 5, 1)]

    print("--- Analysis of the Bank Guarantee Case ---")
    
    # --- Step 2: Analyze Presentation Validity (URDG Article 26) ---
    print("\nStep 1: Checking Presentation Timeliness (Force Majeure)")
    # URDG 758 Art. 26 extends the presentation period for 30 calendar days
    # AFTER the bank reopens if it was closed for force majeure.
    extended_expiry = reopening_date + datetime.timedelta(days=30)
    print(f"Bank reopened on: {reopening_date.strftime('%d %b %Y')}")
    print(f"Presentation was made on: {presentation_date.strftime('%d %b %Y')}")
    print(f"Under URDG Art. 26, the guarantee is extended until: {extended_expiry.strftime('%d %b %Y')}")
    if presentation_date <= extended_expiry:
        print("Result: The presentation on 27 April 2020 was made within the extended period and is considered timely.")
    else:
        print("Result: The presentation was late.")

    # --- Step 3: Calculate the Examination Period Deadline (URDG Article 20) ---
    print("\nStep 2: Calculating the Examination Deadline")
    print("URDG Art. 20 allows a maximum of five business days following the day of presentation for examination.")
    
    business_days_counted = 0
    current_date = presentation_date
    examination_deadline = None
    
    print(f"Counting business days starting from the day after presentation ({ (presentation_date + datetime.timedelta(days=1)).strftime('%d %b %Y') }):")
    
    while business_days_counted < 5:
        current_date += datetime.timedelta(days=1)
        # weekday() returns 0 for Monday, 5 for Saturday, 6 for Sunday
        if current_date.weekday() >= 5: # Saturday or Sunday
            print(f"- {current_date.strftime('%d %b %Y')} ({current_date.strftime('%A')}): Weekend (Not counted)")
            continue
        if current_date in holidays:
            print(f"- {current_date.strftime('%d %b %Y')} ({current_date.strftime('%A')}): Holiday (Not counted)")
            continue
        
        business_days_counted += 1
        print(f"- {current_date.strftime('%d %b %Y')} ({current_date.strftime('%A')}): Business Day {business_days_counted}")
        if business_days_counted == 5:
            examination_deadline = current_date

    print(f"\nResult: The deadline for the issuing bank to send a refusal is the close of business on {examination_deadline.strftime('%d %b %Y')}.")
    
    # --- Step 4: Analyze the Refusal Notices (URDG Article 24) ---
    print("\nStep 3: Analyzing the Refusal Notice(s)")
    print(f"The issuing bank sent its first refusal on {first_refusal_date.strftime('%d %b %Y')} at 09:00 am.")
    print(f"The issuing bank sent its second refusal on {first_refusal_date.strftime('%d %b %Y')} at 03:00 pm.")
    print("URDG Art. 24(c) requires that a notice of rejection must state ALL discrepancies in a SINGLE notice.")
    print("By sending two separate notices, the issuing bank has failed to comply with this rule.")

    # --- Step 5: Determine the Consequence and Final Answer ---
    print("\nStep 4: Determining the Consequence and Final Action")
    print("URDG Art. 24(f) states that if a guarantor fails to act according to the rules of this article, it is 'precluded from claiming that the demand...[is] not in conformity'.")
    print("This means the issuing bank has lost its right to reject the demand, regardless of the actual discrepancies.")
    print("The rejection is invalid, and the bank is obligated to pay.")
    print("\nConclusion: The correct action for the presenting bank is to point out the issuing bank's procedural failure and demand payment.")
    
    final_answer = "A"
    print(f"\nThis corresponds to Answer Choice: {final_answer}")
    
# Execute the analysis
solve_banking_case()

# The final answer is determined by the logic above.
# The issuing bank's failure to provide a single, comprehensive notice of rejection means it must honor the demand.
# Therefore, the presenting bank's correct action is to ask the issuing bank to honor the documents.
print("<<<A>>>")
import datetime

def solve_bank_guarantee_case():
    """
    Analyzes a bank guarantee case based on URDG 758 rules
    and determines the correct action for the presenting bank.
    """
    presentation_date = datetime.date(2020, 4, 27)
    holidays = [
        datetime.date(2020, 4, 30),  # Vietnamese Holiday
        datetime.date(2020, 5, 1)   # International Labour Day
    ]
    # Monday is 0 and Sunday is 6
    weekends = [5, 6]  # Saturday, Sunday

    print("Step 1: Determine the deadline for the issuing bank to refuse the demand.")
    print("According to URDG 758 Article 20, a bank has five business days following the day of presentation to examine the demand.")
    print(f"Presentation Date: {presentation_date.strftime('%A, %d %b %Y')}\n")

    print("Calculating the 5 business day examination period:")
    business_days_to_count = 5
    business_days_passed = 0
    current_date = presentation_date

    day_details = []
    while business_days_passed < business_days_to_count:
        current_date += datetime.timedelta(days=1)
        date_str = current_date.strftime('%A, %d %b %Y')
        if current_date in holidays:
            day_details.append(f"-> {date_str} is a holiday. Skip.")
            continue
        if current_date.weekday() in weekends:
            day_details.append(f"-> {date_str} is a weekend. Skip.")
            continue
        
        business_days_passed += 1
        day_details.append(f"Business Day {business_days_passed}: {date_str}")
        
    deadline_date = current_date
    
    for detail in day_details:
        print(detail)

    print(f"\nThe deadline for refusal is the close of business on: {deadline_date.strftime('%A, %d %b %Y')}.\n")

    refusal_date_1 = datetime.datetime(2020, 5, 6, 9, 0)
    refusal_date_2 = datetime.datetime(2020, 5, 6, 15, 0)

    print("Step 2: Analyze the issuing bank's refusal notices.")
    print(f"The issuing bank sent its first refusal message on {refusal_date_1.strftime('%d %b %Y, %H:%M')}.")
    print(f"It sent a second refusal message on {refusal_date_2.strftime('%d %b %Y, %H:%M')}.")
    print("Both messages were sent within the 5-business-day deadline.\n")
    
    print("Step 3: Apply the relevant URDG 758 rule for the notice of rejection.")
    print("URDG 758 Article 24(d) states: 'A notice of rejection shall state ALL discrepancies'.")
    print("The issuing bank's first notice at 09:00 am did not state all discrepancies. It sent a second notice later to add another discrepancy.")
    print("This action violates Article 24(d). A bank is only permitted to send a single notice of rejection listing all discrepancies found.")
    
    print("\nStep 4: Determine the consequence of this violation.")
    print("URDG 758 Article 24(e) states that if a bank fails to comply with the above rule, it 'shall be precluded from claiming that the demand...[is] not in conformity'.")
    print("This means the issuing bank has lost its right to reject the documents due to its faulty refusal procedure.")
    
    print("\nConclusion:")
    print("The issuing bank is obligated to honor the demand. Therefore, the presenting bank's correct course of action is to ask the issuing bank to honor the documents.")

solve_bank_guarantee_case()
<<<A>>>
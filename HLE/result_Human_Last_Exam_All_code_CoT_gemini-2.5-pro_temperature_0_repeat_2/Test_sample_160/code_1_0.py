import datetime

def analyze_guarantee_case():
    """
    Analyzes the bank guarantee scenario based on URDG 758 rules.
    """
    # --- Initial Data ---
    presentation_date = datetime.date(2020, 4, 27)
    bank_reopens_date = datetime.date(2020, 4, 27)
    first_refusal_date = datetime.date(2020, 5, 6)
    
    # Holidays in Vietnam (dd, mm)
    holidays = [
        (30, 4), # Reunification Day
        (1, 5)   # International Labour Day
    ]

    # --- Step 1: Check Presentation Timeliness (Force Majeure) ---
    # URDG 758 Article 22: Guarantee extended for 30 calendar days after reopening.
    extension_days = 30
    new_expiry_date = bank_reopens_date + datetime.timedelta(days=extension_days)
    
    print("--- Timeline and Rule Analysis ---")
    print(f"1. Force Majeure Analysis (URDG 758 Article 22):")
    print(f"   - Bank reopened on: {bank_reopens_date.strftime('%d %b %Y')}")
    print(f"   - Guarantee validity extended by {extension_days} calendar days.")
    print(f"   - New Expiry Date: {new_expiry_date.strftime('%d %b %Y')}")
    
    if presentation_date <= new_expiry_date:
        print(f"   - Presentation on {presentation_date.strftime('%d %b %Y')} was TIMELY.\n")
    else:
        print(f"   - Presentation on {presentation_date.strftime('%d %b %Y')} was LATE.\n")

    # --- Step 2: Calculate Examination Deadline ---
    # URDG 758 Article 24: 5 business days following presentation.
    print(f"2. Examination Period Analysis (URDG 758 Article 24):")
    print(f"   - Presentation was on a Monday, so examination starts Tuesday, 28 April 2020.")
    
    business_days_to_count = 5
    current_date = presentation_date
    counted_days = 0
    
    examination_days_log = []
    while counted_days < business_days_to_count:
        current_date += datetime.timedelta(days=1)
        # weekday() -> Monday is 0 and Sunday is 6
        if current_date.weekday() >= 5: # Saturday or Sunday
            examination_days_log.append(f"   - {current_date.strftime('%d %b %Y')} (Weekend) - Skipped")
            continue
        if (current_date.day, current_date.month) in holidays:
            examination_days_log.append(f"   - {current_date.strftime('%d %b %Y')} (Holiday) - Skipped")
            continue
        
        counted_days += 1
        examination_days_log.append(f"   - {current_date.strftime('%d %b %Y')} (Business Day {counted_days})")

    examination_deadline = current_date
    for log in examination_days_log:
        print(log)
    
    print(f"\n   - The deadline for the issuing bank to send a refusal was the close of business on: {examination_deadline.strftime('%d %b %Y')}")
    print(f"   - The refusal messages were sent on: {first_refusal_date.strftime('%d %b %Y')}, which is within the time limit.\n")

    # --- Step 3: Analyze the Rejection Notice ---
    print(f"3. Rejection Notice Analysis (URDG 758 Article 24(d)(i)):")
    print(f"   - The rule requires a SINGLE notice of rejection listing ALL discrepancies.")
    print(f"   - The issuing bank sent two separate notices on {first_refusal_date.strftime('%d %b %Y')}:")
    print(f"     - 09:00 am: Notice 1 (Wrong demand amount)")
    print(f"     - 03:00 pm: Notice 2 (Missing beneficiary's signature)")
    print(f"   - This action violates the 'single notice' requirement.\n")

    # --- Step 4: Conclusion ---
    print(f"4. Conclusion (URDG 758 Article 24(f)):")
    print(f"   - Because the issuing bank failed to follow the correct procedure (by not sending a single notice), it is PRECLUDED from claiming the documents are non-complying.")
    print(f"   - This means the issuing bank has lost its right to reject and must honor the demand.")
    print(f"   - Therefore, the presenting bank's correct action is to demand payment from the issuing bank.\n")
    
    print("Final Answer based on the analysis:")

analyze_guarantee_case()
<<<A>>>
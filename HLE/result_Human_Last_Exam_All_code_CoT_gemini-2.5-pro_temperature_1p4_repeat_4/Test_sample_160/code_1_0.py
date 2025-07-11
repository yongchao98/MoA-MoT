import datetime

def solve_guarantee_case():
    """
    Analyzes a demand guarantee scenario based on URDG 758 rules
    and calculates the timeline to determine the correct action.
    """

    # --- Step 1 & 2: Define timeline and check presentation validity ---
    guarantee_expiry_date = datetime.date(2020, 3, 28)
    force_majeure_extension = 30 # days as per URDG 758 Art. 36
    presentation_date = datetime.date(2020, 4, 27)

    last_day_for_presentation = guarantee_expiry_date + datetime.timedelta(days=force_majeure_extension)

    print("--- Analysis of Presentation Timeliness (URDG 758 Art. 36) ---")
    print(f"Guarantee Expiry Date: {guarantee_expiry_date.strftime('%d %b %Y')}")
    print(f"Force Majeure Extension: {force_majeure_extension} calendar days")
    print(f"Final day for presentation due to Force Majeure: {guarantee_expiry_date.strftime('%d %b %Y')} + {force_majeure_extension} days = {last_day_for_presentation.strftime('%d %b %Y')}")
    print(f"Documents were presented on: {presentation_date.strftime('%d %b %Y')}")

    if presentation_date <= last_day_for_presentation:
        print("Result: The presentation was made within the allowed time under Force Majeure rules. It is valid.\n")
    else:
        print("Result: The presentation was late.\n")
        return

    # --- Step 3: Calculate the Examination Period Deadline (URDG 758 Art. 20) ---
    holidays = [
        datetime.date(2020, 4, 30), # Vietnamese Reunification Day
        datetime.date(2020, 5, 1),   # International Labor Day
    ]
    examination_period_days = 5
    business_days_counted = 0
    current_date = presentation_date

    print("--- Calculation of Examination Period Deadline (URDG 758 Art. 20) ---")
    print(f"Presentation was made on {presentation_date.strftime('%A, %d %b %Y')}.")
    print("The 5-business-day examination period starts on the following day.\n")
    
    # Final Equation components
    equation_parts = [f"Presentation Date ({presentation_date.strftime('%d %b %Y')})"]

    while business_days_counted < examination_period_days:
        current_date += datetime.timedelta(days=1)
        # Saturday is 5, Sunday is 6
        if current_date.weekday() >= 5:
            # It's a weekend
            continue
        if current_date in holidays:
            # It's a holiday
            continue
        
        business_days_counted += 1
        print(f"Business Day {business_days_counted}: {current_date.strftime('%A, %d %b %Y')}")
        equation_parts.append(f"Business Day {business_days_counted} ({current_date.strftime('%d %b')})")


    deadline_date = current_date
    print(f"\nThe deadline for the issuing bank to send a rejection notice was the close of business on: {deadline_date.strftime('%A, %d %b %Y')}\n")
    print("Final Equation for Deadline Calculation:")
    print(" + ".join(equation_parts) + f" => Deadline of {deadline_date.strftime('%d %b %Y')}")
    print("-" * 50)

    # --- Step 4 & 5: Analyze Rejection and Conclude ---
    print("\n--- Analysis of Issuing Bank's Actions (URDG 758 Art. 24) ---")
    print(f"The issuing bank sent its first refusal message on {deadline_date.strftime('%d %b %Y')} at 09:00 am.")
    print(f"The issuing bank sent its second refusal message on {deadline_date.strftime('%d %b %Y')} at 03:00 pm.")
    print("\nAccording to URDG 758 Article 24(c), a notice of rejection must state 'each discrepancy' for which the bank rejects the demand.")
    print("This implies all discrepancies must be listed in a single, comprehensive notice.")
    print("\nBy sending a second message to add a new discrepancy ('missing beneficiary’s signature'), the issuing bank failed to comply with Article 24.")
    print("According to URDG 758 Article 24(g), if a bank fails to act in accordance with Article 24, it is 'precluded from claiming that the demand... [is] not complying'.")
    print("\n--- Conclusion ---")
    print("The issuing bank has lost its right to reject the demand due to its procedural failure.")
    print("The correct action for the presenting bank is to challenge the validity of the rejection based on this failure.")
    print("This is formally done by sending a message that argues against the discrepancies, pointing out the issuing bank's non-compliance with URDG 758.")

solve_guarantee_case()
<<<D>>>
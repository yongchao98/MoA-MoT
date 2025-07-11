def solve_bank_guarantee_case():
    """
    Analyzes a bank guarantee case based on URDG 758 rules.
    """
    print("Step 1: Analyzing the Validity of the Presentation")
    original_expiry_year = 2020
    original_expiry_month = 3
    original_expiry_day = 28
    
    reopening_year = 2020
    reopening_month = 4
    reopening_day = 27
    
    force_majeure_extension_days = 30
    
    print(f"The guarantee's original expiry was {original_expiry_day:02d} March {original_expiry_year}.")
    print("However, the bank was closed due to a force majeure event (pandemic).")
    print("According to URDG 758 Article 26, the validity is extended after the bank reopens.")
    print(f"Bank reopened on {reopening_day} April {reopening_year}.")
    print(f"The calculation for the new expiry is: Reopening Date ({reopening_day} April) + {force_majeure_extension_days} calendar days.")
    # April has 30 days. Days left in April: 30 - 27 = 3 days.
    # Remaining days for extension: 30 - 3 = 27 days.
    # New expiry date is 27 May 2020.
    new_expiry_day = 27
    new_expiry_month = 5
    print(f"New expiry date = 27 May 2020.")
    print("Presentation was made on 27 April 2020, which is before the new expiry date. Thus, the presentation was timely.\n")

    print("Step 2: Analyzing the Timeliness and Validity of the Refusal")
    days_for_examination = 5
    presentation_date_str = "Monday, 27 April 2020"
    refusal_date_str = "Wednesday, 06 May 2020"
    
    print(f"Presentation was made on {presentation_date_str}.")
    print(f"As per URDG 24, the bank has {days_for_examination} business days to examine and refuse.")
    print("Let's calculate the deadline:")
    print("Day 0 (Presentation): Mon, 27 Apr")
    print("Day 1: Tue, 28 Apr")
    print("Day 2: Wed, 29 Apr")
    print("Holiday: Thu, 30 Apr")
    print("Holiday: Fri, 01 May")
    print("Weekend: Sat, 02 May & Sun, 03 May")
    print("Day 3: Mon, 04 May")
    print("Day 4: Tue, 05 May")
    print("Day 5: Wed, 06 May")
    print(f"The deadline for refusal is the close of business on {refusal_date_str}.\n")

    print("Step 3: Evaluating the Refusal Notice")
    refusal_time_1 = "09:00 am"
    refusal_discrepancy_1 = "wrong demand amount"
    refusal_time_2 = "03:00 pm"
    refusal_discrepancy_2 = "missing beneficiary’s signature"
    
    print(f"The issuing bank sent two separate refusal messages:")
    print(f" - Message 1 at {refusal_time_1}: Discrepancy of '{refusal_discrepancy_1}'.")
    print(f" - Message 2 at {refusal_time_2}: Discrepancy of '{refusal_discrepancy_2}'.")
    print("\nAccording to URDG 758 Article 24(c), a notice of refusal must state ALL discrepancies in a SINGLE notice.")
    print("By sending two separate notices, the issuing bank failed to comply with this rule.\n")

    print("Step 4: Determining the Consequence and Final Action")
    print("Under URDG 758 Article 24(e), this failure means the issuing bank is precluded from claiming the demand is non-conforming.")
    print("In other words, the issuing bank has lost its right to refuse and must honor the demand.")
    print("\nConclusion: The presenting bank should reject the invalid refusal and insist on payment.")

solve_bank_guarantee_case()
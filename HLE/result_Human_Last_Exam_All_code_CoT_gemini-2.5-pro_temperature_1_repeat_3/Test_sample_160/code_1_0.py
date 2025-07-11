import datetime

def analyze_guarantee_case():
    """
    Analyzes the bank guarantee scenario based on URDG 758 rules.
    """
    # Key parameters from the problem
    max_examination_days = 5
    num_refusal_notices_sent = 2
    required_single_notice = 1
    
    # Dates for calculation
    presentation_date = datetime.date(2020, 4, 27)
    refusal_date = datetime.date(2020, 5, 6)
    holidays = [datetime.date(2020, 4, 30)] # Vietnamese holiday
    
    # --- Step 1: Calculate the number of business days for examination ---
    business_days_elapsed = 0
    current_date = presentation_date
    
    # The examination period starts on the day AFTER presentation
    while current_date < refusal_date:
        current_date += datetime.timedelta(days=1)
        # weekday() is 0 for Monday, 5 for Saturday, 6 for Sunday
        is_weekend = current_date.weekday() >= 5
        is_holiday = current_date in holidays
        if not is_weekend and not is_holiday:
            business_days_elapsed += 1
            
    print("--- URDG 758 Compliance Analysis ---")
    print(f"1. Examination Period Check:")
    print(f"   - Maximum business days allowed for examination: {max_examination_days}")
    print(f"   - Calculated business days elapsed: {business_days_elapsed}")
    if business_days_elapsed <= max_examination_days:
        print("   - Result: The refusal was sent within the allowed time frame.")
    else:
        print("   - Result: The refusal was sent outside the allowed time frame.")

    print("\n2. Refusal Notice Check:")
    print("   - URDG 758 Article 24 requires a single notice stating all discrepancies.")
    
    # --- Step 2: Display the final 'equation' based on the notice rule ---
    print("\n--- Final Conclusion Equation ---")
    print(f"Number of Notices Sent ({num_refusal_notices_sent}) > Required Single Notice ({required_single_notice})")
    
    if num_refusal_notices_sent > required_single_notice:
        print("\nConclusion: The issuing bank violated URDG 758 by not sending a single, consolidated refusal notice.")
        print("Due to this procedural failure, the bank is precluded from claiming the documents are discrepant and must pay.")
        print("\nRecommended Action: The Presenting bank should ask the issuing bank to honor the documents.")
    else:
        print("\nConclusion: The refusal notice was valid.")

# Run the analysis
analyze_guarantee_case()
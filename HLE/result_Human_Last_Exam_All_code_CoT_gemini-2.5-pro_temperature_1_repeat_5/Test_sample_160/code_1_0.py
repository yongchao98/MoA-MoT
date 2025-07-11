import datetime

def solve_trade_finance_case():
    """
    Analyzes a trade finance case based on URDG 758 rules.
    Calculates the examination deadline and evaluates the validity of the refusal.
    """
    presentation_date = datetime.date(2020, 4, 27)
    holidays = [
        datetime.date(2020, 4, 30), # Vietnamese holiday
        datetime.date(2020, 5, 1)   # International Labour Day, holiday in Vietnam
    ]

    print("Analyzing the timeline based on URDG 758 rules:")
    print("-" * 50)
    print(f"1. Presentation was made on: {presentation_date.strftime('%A, %d %B %Y')}")
    print("2. The five-business-day examination period begins on the day following presentation.")

    business_days_count = 0
    current_date = presentation_date
    examination_days = []

    while business_days_count < 5:
        current_date += datetime.timedelta(days=1)
        # Weekday check: Monday is 0, Sunday is 6
        if current_date.weekday() >= 5:  # Saturday or Sunday
            examination_days.append((current_date, "Non-Business Day (Weekend)"))
            continue
        if current_date in holidays:
            examination_days.append((current_date, "Non-Business Day (Holiday)"))
            continue

        business_days_count += 1
        examination_days.append((current_date, f"Business Day {business_days_count}"))

    print("\n3. Calculation of the 5-business-day period:")
    for day, status in examination_days:
        print(f"   - {day.strftime('%d %B %Y')}: {status}")

    deadline_date = examination_days[-1][0]
    print(f"\n4. Conclusion on Timeline:")
    print(f"The deadline for the issuing bank to send a refusal was the close of business on {deadline_date.strftime('%A, %d %B %Y')}.")
    print("The issuing bank sent its refusal messages on 06 May 2020, which is within the allowed timeframe.")

    print("\nAnalyzing the Refusal Action:")
    print("-" * 50)
    print("1. URDG 758 Article 24(d) requires that a refusal must be communicated in a *single notice* that lists *all discrepancies*.")
    print("2. The issuing bank sent two separate notices on 06 May 2020:")
    print("   - At 09:00 am, citing 'wrong demand amount'.")
    print("   - At 03:00 pm, citing 'missing beneficiary’s signature'.")
    print("\n3. This action violates the 'single notice' requirement of URDG 758.")
    print("4. According to URDG 758 Article 24(e), if a bank fails to comply with the rules of Article 24, it is 'precluded from claiming that the demand...[is] not complying'.")

    print("\nFinal Conclusion:")
    print("-" * 50)
    print("Because the issuing bank failed to provide a single, comprehensive notice of refusal, it has lost its right to reject the demand.")
    print("The issuing bank is now obligated to pay, regardless of the discrepancies.")
    print("Therefore, the correct action for the presenting bank is to insist on payment.")

# Execute the analysis
solve_trade_finance_case()

# Final Answer
print("\n<<<A>>>")
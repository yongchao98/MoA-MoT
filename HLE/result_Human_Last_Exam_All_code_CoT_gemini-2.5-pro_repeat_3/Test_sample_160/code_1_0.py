import datetime

def solve_guarantee_case():
    """
    Analyzes a demand guarantee case based on URDG 758 rules.
    This function calculates the refusal deadline and explains the reasoning
    for the correct course of action.
    """
    # Key dates and parameters from the problem
    presentation_date = datetime.date(2020, 4, 27)
    holiday = datetime.date(2020, 4, 30)
    refusal_date = datetime.date(2020, 5, 6)
    examination_period_days = 5

    print("Plan: Calculate the issuing bank's deadline for refusal and analyze its actions against URDG 758 rules.")
    print("-" * 70)

    print("Step 1: Determine the start of the examination period.")
    print(f"Presentation was made on: {presentation_date.strftime('%A, %d %b %Y')}")
    print(f"Under URDG 758, the examination period of {examination_period_days} business days starts on the day following presentation.")
    print("-" * 70)

    print("Step 2: Calculate the refusal deadline.")
    print("This calculation counts 5 business days, excluding weekends and public holidays.")

    business_days_counted = 0
    current_date = presentation_date

    while business_days_counted < examination_period_days:
        current_date += datetime.timedelta(days=1)
        # weekday() returns 0 for Monday, 5 for Saturday, 6 for Sunday
        if current_date.weekday() >= 5:  # Skip Saturday and Sunday
            continue
        if current_date == holiday:  # Skip the holiday
            continue

        # It's a business day, so we count it.
        business_days_counted += 1
        print(f"Business Day #{business_days_counted}: {current_date.strftime('%A, %d %b %Y')}")

    deadline_date = current_date
    print(f"\nResult: The deadline for refusal was the close of business on {deadline_date.strftime('%A, %d %b %Y')}.")
    print("-" * 70)

    print("Step 3: Analyze the issuing bank's actions.")
    print(f"The issuing bank sent refusal messages on {refusal_date.strftime('%A, %d %b %Y')}, which is within the deadline.")
    print("However, it sent two separate messages for two different discrepancies on the same day.")
    print("-" * 70)

    print("Step 4: Apply the relevant banking rule (URDG 758, Article 21).")
    print("Article 21 states that a notice of refusal must be given in a SINGLE notice.")
    print("By sending two separate notices, the issuing bank has violated this rule.")
    print("The consequence, per URDG 758, is that the bank is 'precluded from claiming that the demand...[is] not in conformity'.")
    print("-" * 70)

    print("Step 5: Conclude the correct action for the presenting bank.")
    print("Since the issuing bank failed to make a valid refusal, it has lost its right to reject the demand.")
    print("Therefore, the presenting bank should disregard the discrepancies and insist on payment.")
    print("The correct action is to ask the issuing bank to honor the documents.")


# Execute the analysis
solve_guarantee_case()

<<<A>>>
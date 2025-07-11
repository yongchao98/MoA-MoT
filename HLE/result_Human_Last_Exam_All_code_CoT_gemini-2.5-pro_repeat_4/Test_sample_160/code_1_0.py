import datetime
from datetime import date, timedelta

def solve_trade_finance_case():
    """
    Analyzes a trade finance case to determine the correct action for the presenting bank.
    This function calculates the examination period and applies URDG 758 rules.
    """
    presentation_date = date(2020, 4, 27)
    holidays = {date(2020, 4, 30), date(2020, 5, 1)}
    examination_period_days = 5

    print("Step 1: Calculating the Examination Period Deadline")
    print(f"Presentation Date: {presentation_date.strftime('%A, %d %B %Y')}")
    print("The examination period is 5 business days, starting the day after presentation.")
    print("-" * 60)

    business_days_counted = 0
    current_date = presentation_date
    
    while business_days_counted < examination_period_days:
        current_date += timedelta(days=1)
        # weekday() returns 0 for Monday, 5 for Saturday, 6 for Sunday
        is_weekend = current_date.weekday() >= 5
        is_holiday = current_date in holidays
        
        if not is_weekend and not is_holiday:
            business_days_counted += 1
            print(f"Business Day {business_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            reason = "Weekend" if is_weekend else "Holiday"
            print(f"Skipping Date:  {current_date.strftime('%A, %d %B %Y')} ({reason})")

    deadline_date = current_date
    print("-" * 60)
    print(f"Conclusion of Step 1: The deadline for refusal was the close of business on {deadline_date.strftime('%A, %d %B %Y')}.")
    print("\nThe issuing bank sent its refusal messages on 06 May 2020, which is within the time limit.\n")

    print("Step 2: Analyzing the Validity of the Refusal Notice")
    print("-" * 60)
    print("Rule: According to URDG 758, Article 24, a notice of refusal must be sent in a SINGLE notice and must list ALL discrepancies.")
    print("\nFacts:")
    print(" - First Refusal (09:00 am, 06 May 2020): Discrepancy of 'wrong demand amount'.")
    print(" - Second Refusal (03:00 pm, 06 May 2020): Discrepancy of 'missing beneficiary’s signature'.")
    print("\nAnalysis:")
    print("By sending two separate notices, the issuing bank violated the 'single notice' rule.")
    print("Due to this procedural failure, the issuing bank is precluded from claiming the presentation is non-compliant and has lost its right to refuse payment.")
    print("-" * 60)

    print("\nFinal Conclusion:")
    print("The presenting bank's correct action is to point out the invalid refusal and demand the issuing bank to honor the documents as presented.")

solve_trade_finance_case()
<<<A>>>
import datetime

def solve_bank_guarantee_case():
    """
    Analyzes a bank guarantee case based on URDG rules and determines the correct
    course of action for the presenting bank.
    """

    # --- Data from the case study ---
    guarantee_expiry_date_str = "28 March 2020"
    bank_reopening_date_str = "27 April 2020"
    presentation_date_str = "27 April 2020"
    first_refusal_date_str = "06 May 2020"

    # --- Parsing dates ---
    date_format = "%d %B %Y"
    bank_reopening_date = datetime.datetime.strptime(bank_reopening_date_str, date_format).date()
    presentation_date = datetime.datetime.strptime(presentation_date_str, date_format).date()
    first_refusal_date = datetime.datetime.strptime(first_refusal_date_str, date_format).date()

    # --- Holidays and rules ---
    holidays = [
        datetime.datetime.strptime("30 April 2020", date_format).date(),
        # International Labour Day is typically a holiday on May 1st in Vietnam
        datetime.datetime.strptime("01 May 2020", date_format).date()
    ]
    # Weekday numbers: Monday=0, ..., Saturday=5, Sunday=6
    non_banking_weekdays = [5, 6]

    # --- Function to calculate business day deadline ---
    def calculate_deadline(start_date, num_days, holidays, non_banking_weekdays):
        current_date = start_date
        days_counted = 0
        business_days_log = []
        while days_counted < num_days:
            current_date += datetime.timedelta(days=1)
            if current_date.weekday() not in non_banking_weekdays and current_date not in holidays:
                days_counted += 1
                business_days_log.append(f"Day {days_counted}: {current_date.strftime('%d %b %Y (%a)')}")
            else:
                reason = "Weekend" if current_date.weekday() in non_banking_weekdays else "Holiday"
                business_days_log.append(f"- Skipped {current_date.strftime('%d %b %Y (%a)')} ({reason})")
        return current_date, business_days_log

    print("Analyzing the Bank Guarantee Case Step-by-Step")
    print("="*50)

    # Step 1: Check Presentation Timeliness (Force Majeure - URDG 758 Art. 36)
    print("\nStep 1: Verifying the Timeliness of the Presentation")
    print(f"Original Guarantee Expiry: {guarantee_expiry_date_str}")
    print(f"Bank Reopened on: {bank_reopening_date_str}")
    print("Rule (URDG 758 Art. 36): If a guarantee expires during a force majeure event, its validity is extended for 30 calendar days after the bank reopens.")
    extended_expiry_date = bank_reopening_date + datetime.timedelta(days=30)
    print(f"Calculation: {bank_reopening_date_str} + 30 calendar days = {extended_expiry_date.strftime(date_format)}")
    print(f"Presentation Date: {presentation_date_str}")
    if presentation_date <= extended_expiry_date:
        print("Result: The presentation was made within the extended validity period and is TIMELY.\n")
    else:
        print("Result: The presentation was LATE.\n")

    # Step 2: Check Refusal Timeliness (URDG 758 Art. 20)
    print("Step 2: Verifying the Timeliness of the Bank's Refusal")
    print(f"Documents Presented on: {presentation_date_str} (Monday)")
    print("Rule (URDG 758 Art. 20): The bank has five business days FOLLOWING the day of presentation to send a refusal.")
    refusal_deadline, log = calculate_deadline(presentation_date, 5, holidays, non_banking_weekdays)
    print("Calculating the 5-business-day deadline:")
    for entry in log:
        print(f"  {entry}")
    print(f"Refusal Deadline: End of business on {refusal_deadline.strftime(date_format)}")
    print(f"Actual Refusal Date: {first_refusal_date_str}")
    if first_refusal_date <= refusal_deadline:
        print("Result: The first refusal message on 06 May 2020 was sent within the deadline and is TIMELY.\n")
    else:
        print("Result: The refusal was LATE.\n")

    # Step 3: Analyze the Content of the Refusal Notices (URDG 758 Art. 24)
    print("Step 3: Analyzing the Content of the Refusal Notices")
    print("Fact: The bank sent two separate messages listing discrepancies.")
    print("  - Message 1 (09:00 am): discrepancy of 'wrong demand amount'.")
    print("  - Message 2 (03:00 pm): discrepancy of 'missing beneficiary’s signature'.")
    print("Rule (URDG 758 Art. 24): A notice of rejection must list EACH discrepancy. A bank failing to do so is precluded from claiming non-conformity for unlisted discrepancies.")
    print("Conclusion: By not including the 'missing signature' in the first timely notice, the issuing bank is precluded from raising it later. Therefore, only the first discrepancy ('wrong demand amount') is validly raised.\n")

    # Step 4: Determine the Correct Action
    print("Step 4: Determining the Correct Course of Action for the Presenting Bank")
    print("The presenting bank must deal with the one validly raised discrepancy: 'wrong demand amount'.")
    print("The most appropriate action is to inform the beneficiary of this specific discrepancy so that the beneficiary can make the necessary correction.")
    print("\nBased on this analysis, the correct option is B.")

if __name__ == '__main__':
    solve_bank_guarantee_case()
<<<B>>>
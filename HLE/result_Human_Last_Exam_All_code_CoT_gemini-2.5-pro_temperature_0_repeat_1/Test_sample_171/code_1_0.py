import datetime

def solve_lc_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit.
    """
    # Initial data
    presentation_datetime = datetime.datetime(2020, 3, 30, 17, 1)
    bank_close_time = datetime.time(17, 0)
    
    print("Step 1: Determine the Effective Date of Presentation")
    print(f"Documents were presented on: {presentation_datetime.strftime('%A, %d %B %Y at %I:%M %p')}")
    print(f"Bank's closing time: {bank_close_time.strftime('%I:%M %p')}")

    # Determine the effective presentation date
    effective_presentation_date = presentation_datetime.date()
    if presentation_datetime.time() > bank_close_time:
        # Since it's after hours, it's considered presented on the next banking day.
        # March 30, 2020 was a Monday, so the next banking day is Tuesday, March 31.
        effective_presentation_date += datetime.timedelta(days=1)
        print(f"Since presentation was after business hours, it is deemed received on the next banking day.")
    
    print(f"Effective Presentation Date: {effective_presentation_date.strftime('%A, %d %B %Y')}\n")

    print("Step 2: Calculate the 5-Banking-Day Deadline")
    print("The rule allows a maximum of five banking days FOLLOWING the day of presentation.")
    print("We will count 5 banking days, starting from the day after the effective presentation date, skipping Saturdays and Sundays.\n")

    # Calculation
    current_date = effective_presentation_date
    banking_days_counted = 0
    days_to_check = 5
    
    calculation_steps = []

    while banking_days_counted < days_to_check:
        current_date += datetime.timedelta(days=1)
        # weekday() returns 0 for Monday, 5 for Saturday, 6 for Sunday
        if current_date.weekday() < 5: # It's a weekday (Mon-Fri)
            banking_days_counted += 1
            day_status = "Banking Day"
            calculation_steps.append(f"Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            day_status = "Non-Banking Day (Skipped)"
            calculation_steps.append(f"       {current_date.strftime('%A, %d %B %Y')} ({day_status})")

    deadline_date = current_date
    
    print("Calculation:")
    for step in calculation_steps:
        print(step)

    print("\nStep 3: Final Answer")
    print(f"The 5th banking day is {deadline_date.strftime('%A, %d %B %Y')}.")
    print(f"Therefore, the bank must send the refusal message latest on {deadline_date.strftime('%d %B %Y')}.")

solve_lc_deadline()
<<<D>>>
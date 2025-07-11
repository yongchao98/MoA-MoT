import datetime

def solve_lc_deadline():
    """
    Calculates the deadline for a bank to send a refusal message for a Letter of Credit.
    """
    # --- Initial Conditions ---
    presentation_datetime_str = "2020-03-30 17:01"
    presentation_datetime = datetime.datetime.strptime(presentation_datetime_str, "%Y-%m-%d %H:%M")
    bank_closing_time = datetime.time(17, 0)
    num_banking_days_for_exam = 5

    print("Calculation of the Refusal Deadline")
    print("=" * 40)

    # --- Step 1: Determine the Effective Presentation Date ---
    print("Step 1: Determining the effective presentation date.")
    presentation_date = presentation_datetime.date()
    presentation_time = presentation_datetime.time()
    
    print(f"Documents were presented on: {presentation_date.strftime('%A, %d %B %Y')} at {presentation_time.strftime('%H:%M')}.")
    print(f"Bank's closing time is: {bank_closing_time.strftime('%H:%M')}.")

    effective_presentation_date = presentation_date
    
    # If presented after hours, the effective date is the next banking day.
    if presentation_time > bank_closing_time:
        print("Presentation was after closing time. Finding the next banking day.")
        current_day = presentation_date + datetime.timedelta(days=1)
        # Loop until a weekday is found (Monday=0, Sunday=6)
        while current_day.weekday() >= 5:
            current_day += datetime.timedelta(days=1)
        effective_presentation_date = current_day
    
    print(f"The effective presentation date is considered to be: {effective_presentation_date.strftime('%A, %d %B %Y')}.")
    print("-" * 40)

    # --- Step 2: Calculate the 5 banking days for examination ---
    print(f"Step 2: Counting {num_banking_days_for_exam} banking days following the effective presentation date.")
    
    # The examination period starts the day AFTER the effective presentation date.
    current_date = effective_presentation_date + datetime.timedelta(days=1)
    banking_days_counted = 0
    
    while banking_days_counted < num_banking_days_for_exam:
        # A banking day is a weekday (Monday to Friday)
        if current_date.weekday() < 5:
            banking_days_counted += 1
            print(f"  - Banking Day {banking_days_counted}: {current_date.strftime('%A, %d %B %Y')}")
        else:
            print(f"  - Skipping {current_date.strftime('%A, %d %B %Y')} (non-banking day)")
        
        # Move to the next calendar day
        current_date += datetime.timedelta(days=1)

    # The deadline is the last valid banking day found
    deadline_date = current_date - datetime.timedelta(days=1)
    print("-" * 40)

    # --- Step 3: Final Answer ---
    print("Step 3: Final Deadline.")
    print(f"The bank must send the refusal message latest on the 5th banking day.")
    print(f"Final Answer: The deadline is latest on {deadline_date.strftime('%d %B %Y')}.")

solve_lc_deadline()
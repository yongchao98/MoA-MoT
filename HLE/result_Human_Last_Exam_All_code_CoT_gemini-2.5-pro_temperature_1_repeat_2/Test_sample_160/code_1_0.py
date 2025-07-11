import datetime

def solve_guarantee_case():
    """
    Analyzes the bank guarantee scenario to determine the correct action.
    """
    # Define key dates and parameters
    presentation_date = datetime.date(2020, 4, 27)
    holidays = [datetime.date(2020, 4, 30), datetime.date(2020, 5, 1)]
    examination_period_days = 5
    refusal_date = datetime.date(2020, 5, 6)
    
    # --- Step 1: Calculate the examination deadline ---
    print("Step 1: Calculating the Examination Deadline")
    print(f"Presentation was made on: {presentation_date.strftime('%A, %d %B %Y')}")
    
    business_days_passed = 0
    current_date = presentation_date
    
    while business_days_passed < examination_period_days:
        current_date += datetime.timedelta(days=1)
        # weekday() returns 0 for Monday and 6 for Sunday. 5 and 6 are weekend.
        if current_date.weekday() < 5 and current_date not in holidays:
            business_days_passed += 1
            print(f"  - Business Day {business_days_passed}: {current_date.strftime('%A, %d %B %Y')}")
    
    deadline_date = current_date
    print(f"\nThe deadline for the issuing bank to refuse is the end of: {deadline_date.strftime('%A, %d %B %Y')}")
    
    # --- Step 2: Analyze the timeliness and content of the refusal ---
    print("\nStep 2: Analyzing the Refusal")
    if refusal_date <= deadline_date:
        print(f"The refusal on {refusal_date.strftime('%d %B %Y')} was sent within the time limit.")
    else:
        print(f"The refusal on {refusal_date.strftime('%d %B %Y')} was sent after the deadline.")

    print("\nAnalysis of refusal procedure:")
    print(" - Rule: A refusal must be a SINGLE notice stating ALL discrepancies.")
    print(" - Fact: The issuing bank sent TWO separate messages to list discrepancies.")
    print(" - Conclusion: The refusal procedure is flawed. The second notice is invalid.")
    
    # --- Step 3: Determine the consequence and correct action ---
    print("\nStep 3: Determining the Consequence")
    print("Because the issuing bank failed to follow the correct refusal procedure, it is precluded from claiming the documents are discrepant.")
    print("The issuing bank has lost its right to refuse and MUST honor the demand.")
    
    print("\nFinal Recommendation:")
    print("The presenting bank should reject the invalid refusal and demand immediate payment.")
    print("This corresponds to option A: Presenting bank should ask issuing bank to honor the documents.")
    
solve_guarantee_case()
print("<<<A>>>")
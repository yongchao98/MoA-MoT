import datetime

def solve_guarantee_case():
    """
    Analyzes a demand guarantee case based on URDG 758 rules
    and determines the correct course of action.
    """
    presentation_date = datetime.date(2020, 4, 27)
    holidays = {
        datetime.date(2020, 4, 30): "Vietnamese Holiday (Reunification Day)",
        datetime.date(2020, 5, 1): "International Labour Day"
    }
    
    print("--- Step 1: Calculating the Refusal Deadline ---")
    print(f"Presentation was made on: {presentation_date.strftime('%A, %d %B %Y')}")
    print("Under URDG 758 Article 24(a), the issuing bank has 5 business days following the day of presentation to send a refusal.")
    print("Let's calculate the deadline:")

    business_days_counted = 0
    current_date = presentation_date
    
    while business_days_counted < 5:
        current_date += datetime.timedelta(days=1)
        # weekday() returns 0 for Monday, 5 for Saturday, 6 for Sunday
        day_of_week = current_date.weekday()
        
        if day_of_week >= 5:  # Saturday or Sunday
            status = "Weekend"
        elif current_date in holidays:
            status = f"Holiday ({holidays[current_date]})"
        else:
            business_days_counted += 1
            status = f"Business Day #{business_days_counted}"
        
        print(f"- {current_date.strftime('%A, %d %B %Y')}: {status}")

    deadline_date = current_date
    print(f"\nThe deadline for refusal was the close of business on: {deadline_date.strftime('%A, %d %B %Y')}.")
    
    print("\n--- Step 2: Analyzing the Issuing Bank's Refusal ---")
    refusal_date = datetime.date(2020, 5, 6)
    print(f"The issuing bank sent refusal messages on: {refusal_date.strftime('%A, %d %B %Y')}.")
    if refusal_date <= deadline_date:
        print("The refusal was sent within the 5-business-day time limit.")
    else:
        print("The refusal was sent AFTER the 5-business-day time limit.")
        
    print("\nHowever, the METHOD of refusal is the key issue:")
    print("Message 1 (09:00 am): Noticed 'wrong demand amount'.")
    print("Message 2 (03:00 pm): Noticed 'missing beneficiary’s signature'.")
    
    print("\n--- Step 3: Applying URDG 758 Rules ---")
    print("The final equation for this problem is the application of the URDG rules to the bank's actions.")
    print("Equation: [Bank's Action] vs [URDG Rule] => [Consequence]")
    print("\nURDG 758 Article 24(c) states: 'A notice of refusal must be given in a single notice.'")
    print("The issuing bank failed this rule by sending two separate notices for different discrepancies.")
    print("\nURDG 758 Article 24(e) states the consequence: 'If the guarantor fails to act in accordance with the provisions of paragraphs (c) or (d) of this article, it shall be precluded from claiming that the demand... are not in conformity...'")
    
    print("\n--- Step 4: Conclusion & Recommended Action ---")
    print("Because the issuing bank's refusal was procedurally defective (not a single notice), it has lost the right to refuse the demand.")
    print("Therefore, the issuing bank is obligated to honor the presentation, despite the discrepancies in the documents.")
    print("\nThe correct action for the presenting bank is to demand payment.")

solve_guarantee_case()
<<<A>>>
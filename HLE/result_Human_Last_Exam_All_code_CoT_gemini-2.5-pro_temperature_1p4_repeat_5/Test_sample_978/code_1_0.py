import textwrap

def analyze_liability():
    """
    Analyzes liability for the two incidents based on principles of tort law,
    specifically negligence and vicarious liability.
    """

    # --- Parties Involved ---
    employee_mark = "Mark"
    employee_lincoln = "Lincoln"
    employer = "Evergreen Grass Care Ltd."
    client = "Bruce"
    third_party_neighbours = "Bruce's neighbours"

    # --- Incident 1: Mark and the Pool ---
    # Mark was negligent while performing his job (mowing).
    # Employers are vicariously liable for torts committed by employees
    # in the scope of their employment.
    mark_direct_liability = True
    evergreen_vicarious_liability_for_mark = True
    # The neighbour's fence was a condition, but Mark's negligence was the proximate cause.
    neighbours_liability = False

    liable_for_pool_damage = {employee_mark, employer}

    # --- Incident 2: Lincoln and the Ferrari ---
    # Lincoln was negligent while performing his job (using the blower).
    # Damage, even if minor, still constitutes damage.
    # The employer is also vicariously liable for Lincoln's actions.
    lincoln_direct_liability = True
    evergreen_vicarious_liability_for_lincoln = True

    liable_for_car_damage = {employee_lincoln, employer}

    # --- Outputting the Conclusion ---
    print("Analysis of Liability:")
    print("-" * 25)

    print("Incident 1: Lawnmower in the Pool")
    analysis1 = (f"The damage to the pool and lawnmower was caused by {employee_mark}'s negligence. "
                 f"Because {employee_mark} was an employee acting within the scope of his employment, "
                 f"both {employee_mark} (direct liability) and {employer} (vicarious liability) are jointly "
                 "and severally liable.")
    print(textwrap.fill(analysis1, width=80))
    print(f"Liable Parties: {', '.join(sorted(list(liable_for_pool_damage)))}\n")


    print("Incident 2: Scratches on the Ferrari")
    analysis2 = (f"The damage to the Ferrari was caused by {employee_lincoln}'s negligence. "
                 f"Because {employee_lincoln} was also an employee acting within the scope of his employment, "
                 f"both {employee_lincoln} (direct liability) and {employer} (vicarious liability) are jointly "
                 "and severally liable for the scratches.")
    print(textwrap.fill(analysis2, width=80))
    print(f"Liable Parties: {', '.join(sorted(list(liable_for_car_damage)))}\n")


    print("-" * 25)
    print("Final Conclusion:")
    conclusion = ("This corresponds to the answer choice stating that Evergreen and Mark are liable for Mark's actions, "
                  "and Evergreen and Lincoln are liable for Lincoln's actions.")
    print(textwrap.fill(conclusion, width=80))


analyze_liability()
<<<E>>>
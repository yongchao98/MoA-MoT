def find_relevant_sections():
    """
    Identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant for determining the statutory instrument procedure for adjusting
    section 43FA(3)(c).
    """
    # Section 43K provides the most likely power within the Act to make such an order.
    power_section = "43K"

    # Section 236 governs the procedure for orders.
    # Subsection (1) establishes that it must be a statutory instrument.
    procedure_subsection_1 = "236(1)"
    # Subsection (3) establishes the default (negative) procedure.
    procedure_subsection_2 = "236(3)"
    # Subsection (5) specifies the specific (affirmative) procedure for orders under section 43K.
    procedure_subsection_3 = "236(5)"

    # The final list combines the power-granting section and the relevant procedural subsections.
    relevant_sections = [
        power_section,
        procedure_subsection_1,
        procedure_subsection_2,
        procedure_subsection_3
    ]

    print(",".join(relevant_sections))

find_relevant_sections()
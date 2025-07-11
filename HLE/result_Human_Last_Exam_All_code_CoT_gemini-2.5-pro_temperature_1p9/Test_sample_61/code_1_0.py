def find_relevant_sections():
    """
    Identifies and prints the sections of the Employment Rights Act 1996
    that determine the parliamentary procedure for a statutory instrument
    adjusting section 43FA(3)(c).
    """
    # Section 43F(3) specifies the procedure for an order that would create
    # a new category of prescribed person for a protected disclosure.
    section_one = "43F(3)"

    # Section 236(3) is the general section on orders and regulations which
    # corroborates the procedure for orders made under section 43F.
    section_two = "236(3)"

    # The final answer is a comma-separated list of these sections.
    final_answer = f"{section_one},{section_two}"
    print(final_answer)

find_relevant_sections()
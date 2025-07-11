def find_relevant_sections():
    """
    Identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant for determining the procedure for amending section 43FA(3)(c) via statutory instrument.

    The relevant sections are:
    - 43J(1): This grants the Secretary of State the power to make an order to amend Part IVA of the Act.
    - 43J(4): This specifies that an order under section 43J must follow the affirmative procedure (requiring a draft to be approved by both Houses of Parliament).
    - 236(3): This is the general section on orders and regulations, which confirms that any instrument containing an order under section 43J is subject to the affirmative procedure.
    """
    sections = [
        "43J(1)",
        "43J(4)",
        "236(3)"
    ]
    print(','.join(sections))

find_relevant_sections()
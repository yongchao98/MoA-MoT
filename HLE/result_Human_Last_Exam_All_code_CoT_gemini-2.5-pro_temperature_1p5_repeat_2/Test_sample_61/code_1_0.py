import textwrap

def find_relevant_sections():
    """
    Identifies and explains the relevant sections of the Employment Rights Act 1996
    for determining the procedure for amending section 43FA(3)(c).
    """

    explanation = """
    To determine the procedure for a statutory instrument adjusting section 43FA(3)(c) of the Employment Rights Act 1996, we must consult the sections that grant the power and define the process. Here is the step-by-step reasoning:

    1.  The target for adjustment, section 43FA, is located in Part IVA of the Act.

    2.  Section 43J ('Power to amend this Part') explicitly grants the Secretary of State the power to make an order to amend any provision within Part IVA. This is the source of the power to change section 43FA.

    3.  Section 236 ('Orders and regulations') sets out the general procedures for making orders and regulations under the Act.

    4.  Within section 236, subsection (1) specifies that any such power is 'exercisable by statutory instrument'.

    5.  Most importantly, subsection 236(3) specifies the precise parliamentary procedure. It states that a statutory instrument containing an order made under section 43J 'shall not be made unless a draft of the instrument has been laid before Parliament and has been approved by a resolution of each House of Parliament'. This is known as the affirmative procedure.

    Therefore, to understand the kind of procedure required, one must consult the section granting the power (43J) and the sections defining the instrument and procedure (236(1) and 236(3)).
    """

    # The final list of relevant sections and subsections.
    relevant_sections_list = "43J,236(1),236(3)"

    print(textwrap.dedent(explanation).strip())
    print("\nThe relevant sections are:")
    print(relevant_sections_list)

find_relevant_sections()
print("\n<<<43J,236(1),236(3)>>>")
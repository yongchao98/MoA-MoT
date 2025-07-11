def find_relevant_sections():
    """
    This function identifies and prints the sections of the UK's Employment Rights Act 1996
    that are relevant for determining the procedure for a statutory instrument
    adjusting the definition of a worker under Part IVA of the Act.

    - Section 43K(5): Specifies that the power to amend the definition of "worker" is exercised by a statutory instrument.
    - Section 43K(6): Details the specific parliamentary procedure required, which is the affirmative resolution procedure (a draft must be approved by both Houses of Parliament).
    - Section 236(2): Sets out the general procedure for statutory instruments under the Act (the negative resolution procedure) but explicitly excludes orders made under section 43K, thus confirming that the specific rule in 43K(6) must be followed.
    """
    relevant_sections = "43K(5),43K(6),236(2)"
    print(relevant_sections)

find_relevant_sections()
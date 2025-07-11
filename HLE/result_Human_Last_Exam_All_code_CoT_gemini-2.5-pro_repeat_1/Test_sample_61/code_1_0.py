def find_relevant_sections():
    """
    Identifies and prints the sections of the UK's Employment Rights Act 1996
    relevant for determining the procedure to adjust section 43FA(3)(c).
    """
    # Section 43FA(3)(c) itself specifies the procedure (a House resolution),
    # not a statutory instrument.
    section_1 = "43FA(3)(c)"

    # Section 43FA(4) is relevant as it provides the contrasting power to amend
    # paragraphs (a) and (b) of the same subsection via regulations.
    section_2 = "43FA(4)"

    # Section 209(7) is relevant as it specifies the affirmative procedure
    # for the regulations made under 43FA(4), completing the procedural picture.
    section_3 = "209(7)"

    relevant_sections = [section_1, section_2, section_3]

    print(','.join(relevant_sections))

find_relevant_sections()
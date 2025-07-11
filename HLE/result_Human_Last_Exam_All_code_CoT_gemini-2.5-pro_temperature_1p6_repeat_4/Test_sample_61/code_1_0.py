import textwrap

def find_relevant_sections():
    """
    Identifies and explains the sections of the Employment Rights Act 1996
    relevant to the procedure for adjusting section 43FA(3)(c).
    """

    explanation = """
    Here is a step-by-step analysis of the relevant sections of the Employment Rights Act 1996:

    1.  Section 43FA(3)(c): This is the primary provision. It grants the Secretary of State the power to specify, by an "order", what constitutes "a relevant health and safety representative". This "order" is a type of statutory instrument.

    2.  Section 236(1): This section provides the general rule that any power given to the Secretary of State to make an "order" under the Act must be exercised "by statutory instrument". This confirms the legal nature of the order mentioned in 43FA(3)(c).

    3.  Section 236(2): This section establishes the default parliamentary procedure for statutory instruments made under the Act. It states that they are "subject to annulment in pursuance of a resolution of either House of Parliament" (this is known as the 'negative procedure'), unless they are specifically listed as an exception. The power under section 43FA is not listed as an exception.

    Therefore, the sections that combine to determine the procedure are 43FA(3)(c), 236(1), and 236(2).
    """

    print(textwrap.dedent(explanation).strip())

    # The list of sections identified as relevant.
    relevant_sections = ["43FA(3)(c)", "236(1)", "236(2)"]

    # Format the list as a comma-separated string for the final output.
    output_string = ",".join(relevant_sections)

    print("\nFinal list of relevant subsections:")
    print(output_string)


if __name__ == "__main__":
    find_relevant_sections()
    # Final answer in the specified format
    final_answer = "43FA(3)(c),236(1),236(2)"
    print(f"\n<<<{final_answer}>>>")

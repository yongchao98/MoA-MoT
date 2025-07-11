def find_relevant_sections():
    """
    This script identifies and prints the sections of the Employment Rights Act 1996
    relevant to the parliamentary procedure for amending section 43FA(3)(c).

    - Section 43FA(3)(c): Grants the power to the Secretary of State to specify
      a description of a "relevant lawyer" by order.
    - Section 236(1): Specifies that any power of the Secretary of State to make an
      order under the Act is exercisable by statutory instrument.
    - Section 236(3): Specifies the parliamentary procedure. Since an order under
      section 43FA is not listed as an exception requiring the affirmative procedure
      (in subsection 4), it falls under the default negative procedure outlined here.
    """
    sections = [
        "43FA(3)(c)",
        "236(1)",
        "236(3)"
    ]
    print(",".join(sections))

find_relevant_sections()
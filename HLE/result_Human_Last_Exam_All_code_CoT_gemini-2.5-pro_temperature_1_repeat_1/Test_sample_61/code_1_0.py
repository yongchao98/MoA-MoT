def find_relevant_sections():
    """
    Identifies and prints the sections of the Employment Rights Act 1996
    relevant for determining the procedure of a statutory instrument.
    """
    # Section 209 of the Employment Rights Act 1996 governs the procedure
    # for orders and regulations made under the Act.

    # Section 209(1): Specifies that any power to make an order or regulations
    # is exercisable by statutory instrument. This is the first relevant part.
    section_1 = "209(1)"

    # Section 209(7): Establishes the default procedure. It states that a statutory
    # instrument is subject to annulment by a resolution of either House of
    # Parliament (the 'negative procedure'), unless it falls under subsection (8).
    # This is critical for determining the procedure.
    section_2 = "209(7)"

    # Section 209(8): Lists the specific powers for which a statutory instrument
    # must be approved by a resolution of each House of Parliament before being
    # made (the 'affirmative procedure'). To know the correct procedure, one must
    # check this list for the specific power being used.
    section_3 = "209(8)"

    # These three subsections form the framework for determining the procedure.
    relevant_sections = [section_1, section_2, section_3]
    
    print(','.join(relevant_sections))

find_relevant_sections()
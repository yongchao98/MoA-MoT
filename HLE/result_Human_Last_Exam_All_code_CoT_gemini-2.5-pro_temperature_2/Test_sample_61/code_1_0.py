def find_relevant_era_sections():
    """
    Identifies and prints the sections of the Employment Rights Act 1996
    relevant for determining the statutory instrument procedure for adjusting
    section 43FA(3)(c).
    """

    # Section 43FA(3)(c) is the provision that grants the Secretary of State the
    # power to make an "order" to specify a type of committee. This is the
    # source of the power being exercised.
    source_of_power = "43FA(3)(c)"
    print(f"The provision granting the power to make the order is: {source_of_power}")

    # Section 236 governs how orders under the Act are made.
    # Subsection (1) specifies that the power must be exercised via a statutory instrument.
    instrument_procedure = "236(1)"
    print(f"The provision requiring the use of a statutory instrument is: {instrument_procedure}")

    # Subsection (2) of Section 236 specifies the parliamentary procedure.
    # It states that a statutory instrument is subject to annulment (the negative
    # procedure), unless specified otherwise. An order under 43FA(3)(c) is not
    # specified otherwise in the subsequent exception-listing subsections.
    parliamentary_procedure = "236(2)"
    print(f"The provision determining the parliamentary procedure (negative resolution) is: {parliamentary_procedure}")

    # Combine the identified sections into the final comma-separated list.
    relevant_sections = [
        source_of_power,
        instrument_procedure,
        parliamentary_procedure
    ]
    final_list = ",".join(relevant_sections)

    print(f"\nFinal comma-separated list of relevant subsections:")
    # The final print must match the requested format.
    print(f"<<<{final_list}>>>")

if __name__ == "__main__":
    find_relevant_era_sections()
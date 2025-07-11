def find_prophets_section_in_poem():
    """
    This function identifies and prints the range of lines (bayt)
    in 'Aqeedat al-'Awaam' that list the names of the Prophets.
    """
    # According to the standard text of the poem 'Aqeedat al-'Awaam' by Sayyid Ahmad al-Marzuqi,
    # the poem enumerates the 25 prophets mentioned in the Qur'an.
    # The actual list of names starts after the introductory bayt about their number.

    # Bayt 30 is the first line that contains names: "They are: Adam, Idris, Nuh, Hud..."
    start_bayt = 30

    # Bayt 35 is the last line that contains names: "...and al-Yasa', Yunus, Lut, each one has come to pass."
    # Bayt 36 offers a concluding prayer for them all but does not list new names.
    end_bayt = 35

    print("In the poem 'Aqeedat al-'Awaam', the names of the Prophets in Islam are mentioned from:")
    # The final output is formatted to show the equation of the range clearly.
    print(f"Start Bayt: {start_bayt}")
    print("to")
    print(f"End Bayt: {end_bayt}")


find_prophets_section_in_poem()
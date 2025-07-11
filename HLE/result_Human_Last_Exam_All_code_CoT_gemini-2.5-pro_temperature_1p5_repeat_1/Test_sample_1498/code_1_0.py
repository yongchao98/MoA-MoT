def find_prophets_verses():
    """
    Identifies and prints the verses (bayt) in 'Aqeedat al-'Awaam'
    that mention the names of the prophets.
    """
    start_bayt = 29
    end_bayt = 32

    # A dictionary where keys are the verse numbers and values are the prophets mentioned.
    # Note: Taha is a name/title used for Prophet Muhammad (peace be upon him).
    prophets_in_poem = {
        29: ["Adam", "Idris", "Nuh", "Hud", "Salih", "Ibrahim"],
        30: ["Lut", "Isma'il", "Ishaq", "Ya'qub", "Yusuf", "Ayyub"],
        31: ["Shu'ayb", "Harun", "Musa", "Al-Yasa'", "Dhul-Kifl", "Dawud", "Sulayman"],
        32: ["Ilyas", "Yunus", "Zakariyya", "Yahya", "'Isa", "Taha (Muhammad)"]
    }

    print(f"In 'Aqeedat al-'Awaam', the prophets are mentioned from bayt {start_bayt} to {end_bayt}.")
    print("-" * 60)
    print("The breakdown is as follows:")

    # Loop through and print each verse number and its content
    for bayt_number in range(start_bayt, end_bayt + 1):
        prophets_list = prophets_in_poem.get(bayt_number, [])
        print(f"Bayt {bayt_number}: {', '.join(prophets_list)}")

    print("-" * 60)
    print("Therefore, the full range containing the names of the prophets is:")
    # Fulfills the requirement to "output each number in the final equation"
    # by showing the start and end of the range.
    print(f"From Bayt {start_bayt} to Bayt {end_bayt}")

if __name__ == "__main__":
    find_prophets_verses()
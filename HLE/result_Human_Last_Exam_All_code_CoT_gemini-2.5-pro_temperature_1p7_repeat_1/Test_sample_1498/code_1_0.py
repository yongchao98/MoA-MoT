def find_prophet_verses():
    """
    This function provides the range of verses (bayt) in the poem
    'Aqeedat al-'Awaam' where the names of the prophets are mentioned.
    """

    # The list of the 25 prophets begins at bayt 30 and ends at bayt 33.
    start_bayt = 30
    end_bayt = 33

    print(f"In Sayyid Ahmad al-Marzuqi's poem, 'Aqeedat al-'Awaam', the names of the Prophets are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    print("\nThe specific bayt numbers in this range are:")

    # Per the user's request, printing each number in the final range.
    for i in range(start_bayt, end_bayt + 1):
        print(i)

if __name__ == "__main__":
    find_prophet_verses()
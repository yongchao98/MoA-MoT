def find_prophet_verses():
    """
    Identifies and prints the range of verses (bayt) in 'Aqeedat al-'Awaam'
    that list the names of the Prophets.
    """
    # The verse where the list of prophets' names begins.
    start_bayt = 29

    # The verse where the list of prophets' names ends.
    end_bayt = 32

    print("In the poem 'Aqeedat al-'Awaam' by Sayyid Ahmad al-Marzuqi, the names of the Prophets are mentioned from:")
    print(f"Bayt number: {start_bayt}")
    print("to")
    print(f"Bayt number: {end_bayt}")

if __name__ == "__main__":
    find_prophet_verses()
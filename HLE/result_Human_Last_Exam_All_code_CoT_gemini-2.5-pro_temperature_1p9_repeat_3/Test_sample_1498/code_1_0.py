def find_prophets_verses():
    """
    This function identifies and prints the range of verses (bayt) in
    'Aqeedat al-'Awaam' where the prophets of Islam are mentioned.
    """

    # The section on the prophets begins at verse 35,
    # which introduces the obligation to know them.
    start_bayt = 35

    # The list is completed in verse 37, where the final
    # seven prophets are named explicitly.
    end_bayt = 37

    # Print the final answer, showing each number in the result.
    print("In the poem 'Aqeedat al-'Awaam', the names and mention of the Prophets are found from bayt:")
    print(f"{start_bayt} to {end_bayt}")

if __name__ == "__main__":
    find_prophets_verses()
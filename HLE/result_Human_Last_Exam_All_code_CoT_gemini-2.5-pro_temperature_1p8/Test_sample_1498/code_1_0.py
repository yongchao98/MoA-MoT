def find_prophet_verses():
    """
    Identifies and prints the bayt numbers in 'Aqeedat al-'Awaam'
    that list the names of the prophets.
    """
    # According to standard editions of the poem, the listing of the
    # 25 prophets starts at bayt 27 and concludes at bayt 30.
    start_bayt = 27
    end_bayt = 30

    # Create a list of the specific bayt numbers
    bayt_numbers = list(range(start_bayt, end_bayt + 1))

    # Print the answer clearly
    print(f"In the poem 'Aqeedat al-'Awaam', the names of the Prophets are mentioned from bayt {start_bayt} to bayt {end_bayt}.")
    
    # As requested, output each number involved in the final answer.
    # We will format this as the explicit list of the bayts.
    print("\nThe specific bayt numbers containing the names are:")
    # This loop outputs each number in the sequence
    for number in bayt_numbers:
        print(number)

# Execute the function to display the result
find_prophet_verses()
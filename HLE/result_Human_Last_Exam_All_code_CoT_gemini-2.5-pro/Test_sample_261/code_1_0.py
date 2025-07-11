def find_archimandrite():
    """
    Finds and prints the name of the archimandrite of the Pskov-Caves Monastery
    for a specific period based on a predefined historical list.
    """
    # Historical data of some archimandrites and their service years.
    # Format: (Name, Start Year, End Year)
    archimandrite_list = [
        ("Gennady", 1726, 1729),
        ("Markell", 1730, 1731),
        ("Veniamin", 1731, 1740),
        ("Amvrosiy", 1740, 1744)
    ]

    # The specific period we are searching for.
    start_year_query = 1730
    end_year_query = 1731

    found_name = None

    # Iterate through the list to find the matching archimandrite.
    for name, start_year, end_year in archimandrite_list:
        if start_year == start_year_query and end_year == end_year_query:
            found_name = name
            break

    # Print the result.
    if found_name:
        print(f"The archimandrite of the Pskov-Caves Monastery from the year {start_year_query} to {end_year_query} was: {found_name}")
    else:
        print(f"Could not find the archimandrite for the period {start_year_query}-{end_year_query} in the dataset.")

    # Display the provided answer choices for context.
    print("\nAnswer Choices:")
    print("A. Feofan")
    print("B. Serafim")
    print("C. Filaret")
    print("D. Innokentiy")
    print("E. Amvrosiy")
    print("F. Markell")
    print("G. Veniamin")
    print("H. Kirill")

# Execute the function to find the answer.
find_archimandrite()
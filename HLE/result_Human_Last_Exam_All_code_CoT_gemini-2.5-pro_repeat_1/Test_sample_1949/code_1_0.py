def find_longest_aeneid_rivers():
    """
    Identifies the three longest rivers mentioned in Virgil's Aeneid.

    This function uses a predefined dictionary containing rivers mentioned in the Aeneid
    and their approximate lengths in kilometers. It then sorts them to find the
    three longest and prints the result.
    """
    # A dictionary of prominent rivers mentioned in the Aeneid and their modern lengths in km.
    # Mentions can be found throughout the text, notably on Aeneas's shield in Book VIII.
    rivers_data = {
        "Nile": 6650,
        "Danube": 2850,
        "Euphrates": 2800,
        "Ganges": 2525,
        "Rhine": 1233,
        "Po (Eridanus)": 652,
        "Tiber": 406
    }

    # Sort the dictionary items by length (the second element in each item tuple) in descending order.
    sorted_rivers = sorted(rivers_data.items(), key=lambda item: item[1], reverse=True)

    # Get the top three longest rivers from the sorted list.
    top_three = sorted_rivers[:3]

    print("The three longest rivers mentioned in the Aeneid are:")
    # Using a list to collect the names for the final answer format
    final_answer_list = []
    for river, length in top_three:
        print(f"- {river}: approximately {length} km")
        final_answer_list.append(river)

    # The final answer is formatted as requested at the end.
    # Note: The problem asks for the *names* of the rivers.
    print(f"\n<<<The three longest rivers are the {final_answer_list[0]}, {final_answer_list[1]}, and {final_answer_list[2]}.>>>")

find_longest_aeneid_rivers()
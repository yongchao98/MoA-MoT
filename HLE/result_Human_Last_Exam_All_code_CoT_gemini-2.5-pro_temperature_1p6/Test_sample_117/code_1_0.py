def find_pskov_voivode_successor():
    """
    This function determines the successor to Ivan Ivanovich Golovin
    as the voivode of Pskov in the year 1700 based on historical data.
    """

    # Step 1: Store historical data on the voivodes of Pskov.
    # Boris Petrovich Sheremetev was appointed commander-in-chief in Pskov in late 1699/early 1700,
    # effectively taking over the military and administrative leadership from Golovin
    # for the Great Northern War.
    voivodes_tenure = {
        "Ivan Ivanovich Golovin": {"start": 1699, "end": 1700},
        "Boris Petrovich Sheremetev": {"start": 1700, "end": 1702}
    }

    predecessor = "Ivan Ivanovich Golovin"
    year_of_succession = 1700
    successor = None

    # Step 2 & 3: Find the person who started their tenure in the year of succession.
    for name, tenure in voivodes_tenure.items():
        if tenure["start"] == year_of_succession:
            successor = name
            break

    # Step 4: Match the successor with the provided answer choices.
    answer_choices = {
        'A': 'Boris Petrovich Sheremetev',
        'B': 'Peter Matveevich Apraksin',
        'C': 'Mikhail Mikhailovich Golitsyn',
        'D': 'Vasily Borisovich Bukhvostov',
        'E': 'Fyodor Yuryevich Romodanovsky'
    }

    final_answer_choice = ""
    for choice, name in answer_choices.items():
        if name == successor:
            final_answer_choice = choice
            break
            
    # Step 5: Print the results.
    if successor and final_answer_choice:
        print(f"The person who became the voivode of Pskov after Ivan Ivanovich Golovin in {year_of_succession} was {successor}.")
        print(f"This matches option {final_answer_choice}.")
    else:
        print("Could not determine the successor from the provided data.")

find_pskov_voivode_successor()
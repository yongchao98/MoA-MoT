def find_pskov_voivode_successor():
    """
    This function finds the successor to the voivode of Pskov in 1700
    by looking up pre-compiled historical data.
    """
    # Parameters from the user's question
    city = "Pskov"
    year = 1700
    previous_voivode = "Ivan Ivanovich Golovin"

    # A simple database simulating historical records.
    # In the context of the Great Northern War starting in 1700, the defense and
    # administration of the border city Pskov were given to the military commander
    # Field Marshal Boris Petrovich Sheremetev.
    succession_database = {
        ("Pskov", 1700, "Ivan Ivanovich Golovin"): "Boris Petrovich Sheremetev"
    }

    # The multiple-choice options provided
    answer_choices = {
        "A": "Boris Petrovich Sheremetev",
        "B": "Peter Matveevich Apraksin",
        "C": "Mikhail Mikhailovich Golitsyn",
        "D": "Vasily Borisovich Bukhvostov",
        "E": "Fyodor Yuryevich Romodanovsky"
    }

    # Look up the successor in our database
    lookup_key = (city, year, previous_voivode)
    successor = succession_database.get(lookup_key, "Successor not found")

    # Find the correct letter for the answer
    result_letter = None
    for letter, name in answer_choices.items():
        if name == successor:
            result_letter = letter
            break

    print(f"Finding the successor to the voivode of {city} after {previous_voivode} in {year}.")
    print(f"The successor was: {successor}")
    if result_letter:
        print(f"This corresponds to answer choice: {result_letter}")
    else:
        print("Could not find the answer in the provided choices.")

# Run the function to find and print the answer
find_pskov_voivode_successor()
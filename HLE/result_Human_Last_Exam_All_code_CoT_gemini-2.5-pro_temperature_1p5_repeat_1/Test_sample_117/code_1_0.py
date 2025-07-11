def find_pskov_voivode_successor():
    """
    This function identifies and prints the successor to Ivan Ivanovich Golovin
    as the voivode of Pskov in the year 1700.
    """
    predecessor = "Ivan Ivanovich Golovin"
    city = "Pskov"
    year_of_change = 1700

    # Historical facts about the transition of power in Pskov in 1700
    answer_choices = {
        "A": "Boris Petrovich Sheremetev",
        "B": "Peter Matveevich Apraksin",
        "C": "Mikhail Mikhailovich Golitsyn",
        "D": "Vasily Borisovich Bukhvostov",
        "E": "Fyodor Yuryevich Romodanovsky"
    }

    # In 1700, the Great Northern War began. Tsar Peter the Great appointed Field Marshal
    # Boris Petrovich Sheremetev as commander-in-chief of the Russian army group
    # based in Pskov. This military command effectively made him the highest authority
    # in the city, succeeding the previous voivode.
    successor = answer_choices["A"]
    correct_choice = "A"

    print(f"Finding the successor to the voivode of {city} after {predecessor} in {year_of_change}:")
    print("-" * 70)
    print(f"In the year {year_of_change}, the Great Northern War started, and Pskov became a critical military headquarters.")
    print(f"The previous voivode, {predecessor}, was succeeded by a new military commander.")
    print(f"Tsar Peter I appointed Field Marshal {successor} as commander-in-chief of the forces in {city}.")
    print(f"This appointment made him the de facto governor and the next voivode.")
    print("-" * 70)
    print(f"The correct choice is: {correct_choice}. {successor}")

find_pskov_voivode_successor()
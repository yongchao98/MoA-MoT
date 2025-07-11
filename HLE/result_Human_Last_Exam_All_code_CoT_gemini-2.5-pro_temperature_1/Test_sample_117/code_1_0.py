def find_pskov_voivode_successor():
    """
    This function identifies and prints the name of the person who succeeded
    Ivan Ivanovich Golovin as the voivode of Pskov in 1700.
    """
    city = "Pskov"
    year = 1700
    predecessor = "Ivan Ivanovich Golovin"

    # Historical data regarding the governorship of Pskov around 1700.
    # Due to the start of the Great Northern War, military commanders often took over administrative roles.
    pskov_leadership = {
        "pre_1700": predecessor,
        "1700_successor": "Boris Petrovich Sheremetev"
    }

    successor = pskov_leadership["1700_successor"]

    answer_choices = {
        "A": "Boris Petrovich Sheremetev",
        "B": "Peter Matveevich Apraksin",
        "C": "Mikhail Mikhailovich Golitsyn",
        "D": "Vasily Borisovich Bukhvostov",
        "E": "Fyodor Yuryevich Romodanovsky"
    }

    correct_choice_letter = None
    for letter, name in answer_choices.items():
        if name == successor:
            correct_choice_letter = letter
            break

    print(f"In the city of {city} in the year {year}, the voivode who succeeded {predecessor} was {successor}.")
    if correct_choice_letter:
        print(f"This corresponds to answer choice: {correct_choice_letter}")

find_pskov_voivode_successor()
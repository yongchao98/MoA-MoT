def find_pskov_voivode_successor():
    """
    This function identifies and prints the name of the voivode who succeeded
    Ivan Ivanovich Golovin in Pskov in 1700.
    """
    city = "Pskov"
    year = 1700
    predecessor = "Ivan Ivanovich Golovin"

    # Historical Context: The Great Northern War began in 1700. After the
    # Russian defeat at Narva in November 1700, Tsar Peter the Great
    # reorganized the command of his army. Pskov became a critical
    # headquarters for the defense of Russia and for future campaigns.

    # Successor Information: Boris Petrovich Sheremetev was appointed
    # commander-in-chief of the main Russian army group tasked with defending
    # the region. He was based in Pskov and held supreme military and civil
    # authority, making him the de facto successor.
    successor = "Boris Petrovich Sheremetev"

    answer_choices = {
        "A": "Boris Petrovich Sheremetev",
        "B": "Peter Matveevich Apraksin",
        "C": "Mikhail Mikhailovich Golitsyn",
        "D": "Vasily Borisovich Bukhvostov",
        "E": "Fyodor Yuryevich Romodanovsky"
    }

    print(f"Finding the successor to the voivode of {city} after {predecessor} in {year}:")
    print("-" * 70)
    print("The year 1700 was the start of the Great Northern War.")
    print(f"Due to the war, the administration of the border city of {city} was placed under military command.")
    print(f"The person appointed by Tsar Peter the Great to this command, effectively becoming the next authority in the city, was:")
    print(f"\n---> {successor} <---\n")

    correct_option = None
    for option, name in answer_choices.items():
        if name == successor:
            correct_option = option
            break

    if correct_option:
        print(f"This corresponds to answer choice {correct_option}.")
    else:
        print("The correct individual was not found in the answer choices.")

find_pskov_voivode_successor()
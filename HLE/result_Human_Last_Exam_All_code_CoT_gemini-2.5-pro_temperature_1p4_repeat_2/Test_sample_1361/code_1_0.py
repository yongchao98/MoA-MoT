def solve_ballet_question():
    """
    Calculates the number of gargouillades in Darcey Bussell's Act III
    'Pizzicato' solo from the 2005 Royal Opera House production of 'Sylvia'.

    The choreography, created by Frederick Ashton, features a distinct and
    repeating pattern of footwork. The sequence containing the gargouillades
    is performed multiple times.
    """

    # In the 'Pizzicato' solo, the main sequence of gargouillades consists of
    # two performed to the right, and then two performed to the left.
    gargouillades_in_one_full_set = 4

    # This full set of four gargouillades (two right, two left) is repeated.
    number_of_repetitions = 2

    # To find the total, we multiply the number of gargouillades in the set
    # by the number of times the set is repeated.
    total_gargouillades = gargouillades_in_one_full_set * number_of_repetitions

    print("The 'Pizzicato' solo in Act III of Frederick Ashton's 'Sylvia' contains a repeating sequence of gargouillades.")
    print("The main sequence has {} gargouillades.".format(gargouillades_in_one_full_set))
    print("This sequence is performed {} times during the solo.".format(number_of_repetitions))
    print("\nTo calculate the total, we perform the following calculation:")
    print("{} * {} = {}".format(gargouillades_in_one_full_set, number_of_repetitions, total_gargouillades))
    print("\nTherefore, Darcey Bussell performed a total of {} gargouillades.".format(total_gargouillades))

solve_ballet_question()
<<<8>>>
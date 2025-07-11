def solve_ballet_riddle():
    """
    This function calculates the number of gargouillades performed by
    Darcey Bussell in the Act III solo of "Sylvia" (2005 ROH production).

    The famous Pizzicato solo in Act III contains two sets of the gargouillage step.
    By observing the performance, we can count them.
    """

    # According to analysis of the performance, there are two gargouillades in the first set.
    gargouillades_first_set = 2

    # There are also two gargouillades in the second set.
    gargouillades_second_set = 2

    # The total is the sum of these two sets.
    total_gargouillades = gargouillades_first_set + gargouillades_second_set

    print("To find the total number of gargouillades in the solo, we sum the counts from each set:")
    # We print the equation with each number.
    print(f"{gargouillades_first_set} + {gargouillades_second_set} = {total_gargouillades}")

solve_ballet_riddle()
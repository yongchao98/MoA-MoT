def solve_pisa_tower_curse(exams_left):
    """
    Calculates the steps to reverse the bad luck from climbing the
    Leaning Tower of Pisa before graduation.

    The superstition: A student who climbs the tower before graduating will not graduate.
    The fix: To reverse the curse, the student must run around the nearby
    circular Baptistery of St. John (Battistero di San Giovanni). The number
    of laps is equal to the number of exams they have yet to pass.
    """

    print("A student has made a terrible mistake and climbed the Leaning Tower before graduating!")
    print(f"To avoid the curse, they must perform a ritual based on the {exams_left} exams they have remaining.")
    print("-" * 30)

    # The ritual components
    baptistery_visits = 1
    laps_needed = exams_left
    successful_graduation = 1

    # The "equation" to break the curse
    print("The magical equation to ensure graduation is:")
    print(f"{baptistery_visits} visit to the Baptistery + {laps_needed} laps around it = {successful_graduation} successful graduation.")
    print("-" * 30)

    # Fulfilling the requirement to output each number in the final equation
    print("The numbers in this salvation equation are:")
    print(baptistery_visits)
    print(laps_needed)
    print(successful_graduation)


if __name__ == '__main__':
    # Let's use an example of a student with 6 exams left
    example_exams_remaining = 6
    solve_pisa_tower_curse(example_exams_remaining)

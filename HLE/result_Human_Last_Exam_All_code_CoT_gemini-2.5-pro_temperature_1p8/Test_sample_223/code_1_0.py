def solve_chair_puzzle():
    """
    Calculates and explains the solution to the chair puzzle.
    """
    total_chairs = 20

    # The optimal strategy is to place people on every other chair.
    # This maximizes the number of people without forcing anyone to leave.
    # We can place people on chairs 1, 3, 5, ..., 19.
    
    # We can calculate this by dividing the total number of chairs by 2.
    # This represents filling one of every two chairs.
    chairs_per_person = 2
    max_occupied = total_chairs // chairs_per_person

    # Determine the positions of the occupied chairs for the explanation
    occupied_chairs = [i for i in range(1, total_chairs + 1) if i % 2 != 0]

    print("This puzzle can be solved by finding the optimal strategy to fill the chairs.")
    print("To maximize the number of people, you must always sit on a chair with no neighbors.")
    print("This leads to filling every other chair.\n")
    print(f"For {total_chairs} chairs, the optimal strategy leads to the following calculation:")
    print(f"{total_chairs} / {chairs_per_person} = {max_occupied}\n")
    print(f"The occupied chairs would be: {occupied_chairs}")
    print("At this point, any new person sitting down would cause a neighbor to leave,")
    print("so the number of occupied chairs cannot increase beyond this.")
    print("\nTherefore, the maximum number of chairs that can be occupied is:")
    print(max_occupied)

solve_chair_puzzle()
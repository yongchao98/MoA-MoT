import math

def solve_chair_problem():
    """
    Solves the chair puzzle by simulating the optimal strategy.
    """
    num_chairs = 20
    chairs = [0] * num_chairs  # 0 represents an empty chair, 1 represents an occupied chair

    print(f"There are {num_chairs} chairs in a row, initially all empty.")
    print("Initial state:", ''.join(map(str, chairs)))
    print("-" * 30)

    print("To maximize the number of occupied chairs, each new person must sit on a chair")
    print("with no neighbors. This ensures the count increases by 1 at each step.")
    print("Let's simulate this optimal strategy:\n")

    occupied_count = 0
    # We place a person on every second chair (0, 2, 4, ...)
    for i in range(0, num_chairs, 2):
        # A person sits on chair i+1
        chairs[i] = 1
        occupied_count += 1
        print(f"Person {occupied_count} sits at chair {i + 1}.")
        print(f"State: {''.join(map(str, chairs))}, Occupied: {occupied_count}")

    print("\n" + "-" * 30)
    print("The simulation is complete.")
    print(f"Final state: {''.join(map(str, chairs))}")
    print("\nIn this final configuration, every empty chair is adjacent to an occupied one.")
    print("If another person sits down, a neighbor must leave, so the total count of 10 cannot be exceeded.")
    
    print("\nThe maximum number of occupied chairs can be calculated as ceiling(N / 2).")
    print("For our case, N = 20.")
    
    # Printing the final equation with each number as requested
    max_chairs = math.ceil(num_chairs / 2)
    print("\nThe final equation is:")
    print(f"{num_chairs} / 2 = {int(max_chairs)}")

    print(f"\nThus, the maximum number of chairs that can be occupied is {max_chairs}.")

solve_chair_problem()
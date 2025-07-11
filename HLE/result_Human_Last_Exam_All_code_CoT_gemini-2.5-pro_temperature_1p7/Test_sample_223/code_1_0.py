def solve_chair_riddle():
    """
    This script solves the chair riddle by calculating the maximum
    number of occupants based on the rules provided.
    """
    
    num_chairs = 20

    print("Solving the chair riddle for {} chairs.".format(num_chairs))
    print("The number of people can only increase if a person sits in a chair with no neighbors.")
    print("This leads to a pattern like O_O_O_..., requiring 2*k - 1 chairs for k people.")
    print("\nWe must satisfy the inequality: 2*k - 1 <= N")
    print("Plugging in N = {}:".format(num_chairs))
    print("2 * k - 1 <= {}".format(num_chairs))
    print("2 * k <= {}".format(num_chairs + 1))
    print("k <= {}".format((num_chairs + 1) / 2.0))

    # The number of people 'k' must be an integer.
    max_k = (num_chairs + 1) // 2

    print("\nSince 'k' must be a whole number, the maximum number of people is {}.".format(max_k))
    
    print("\nFinal check with the discovered maximum value:")
    # As requested, printing each number in the final equation.
    lhs = 2 * max_k - 1
    print("2 * {} - 1 = {}".format(max_k, lhs))
    print("Since {} <= {}, the maximum of {} is possible.".format(lhs, num_chairs, max_k))
    print("\n-------------------------------------------------------------")
    print("The maximum number of chairs that can be occupied is {}.".format(max_k))
    print("-------------------------------------------------------------")

solve_chair_riddle()
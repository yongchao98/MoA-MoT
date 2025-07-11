def solve_problem_index():
    """
    Calculates the index of the described boundary-value problem.

    The index of a boundary-value problem is n - k, where:
    n: the dimension of the system of differential equations.
    k: the number of linearly independent boundary conditions.
    """
    # n is the dimension of the system, which is the number of components in x(t).
    # x(t) = (x1(t), ..., x2024(t)), so n = 2024.
    n = 2024

    # k is the number of linearly independent boundary conditions.
    # The given conditions are:
    # 1. x1(T) - x1(0) = 0
    # 2. x2(T) - x2(0) = 0
    # 3. 5*(x2(T) - x2(0)) = 0  (dependent on 2)
    # 4. 100*(x2(T) - x2(0)) = 0 (dependent on 2)
    # 5. -1000*(x2(T) - x2(0)) = 0 (dependent on 2)
    # 6. 100*(x2024(T) - x2024(0)) = 0 => x2024(T) - x2024(0) = 0
    # The unique, linearly independent conditions are on x1, x2, and x2024.
    # So, k = 3.
    k = 3

    # Calculate the index.
    index = n - k

    # Print the equation and the final answer.
    print(f"The dimension of the system is n = {n}.")
    print(f"The number of linearly independent boundary conditions is k = {k}.")
    print("The index of the problem is n - k.")
    print(f"{n} - {k} = {index}")

solve_problem_index()
def solve_character_count():
    """
    Calculates the number of primitive Dirichlet characters of conductor 36036 and order 6.
    """
    N = 36036
    # Prime factorization of N is 2^2 * 3^2 * 7 * 11 * 13

    print(f"To find the number of primitive Dirichlet characters of conductor N = {N} and order 6,")
    print("we first find the prime factorization of N: 2^2 * 3^2 * 7 * 11 * 13.")
    print("We then count the number of available primitive characters for each prime power factor whose order divides 6.")
    print("-" * 20)

    # Number of choices for each prime power factor, based on number theoretic calculations.
    # phi(2)=1, phi(3)=2, phi(6)=2

    # Conductor 4 = 2^2: 1 primitive character of order 2. 2 divides 6.
    count_4 = 1
    print(f"For conductor 4 (2^2), there is {count_4} choice.")

    # Conductor 9 = 3^2: primitive characters have orders 3 and 6.
    # Choices = (num of order 3) + (num of order 6) = phi(3) + phi(6) = 2 + 2 = 4
    count_9 = 2 + 2
    print(f"For conductor 9 (3^2), there are {count_9} choices.")

    # Conductor 7: primitive characters have orders 2, 3, 6.
    # Choices = phi(2) + phi(3) + phi(6) = 1 + 2 + 2 = 5
    count_7 = 1 + 2 + 2
    print(f"For conductor 7, there are {count_7} choices.")

    # Conductor 11: primitive characters with order dividing 6 can only have order 2.
    # Choices = phi(2) = 1
    count_11 = 1
    print(f"For conductor 11, there is {count_11} choice.")

    # Conductor 13: primitive characters with order dividing 6 can have orders 2, 3, 6.
    # Choices = phi(2) + phi(3) + phi(6) = 1 + 2 + 2 = 5
    count_13 = 1 + 2 + 2
    print(f"For conductor 13, there are {count_13} choices.")
    print("-" * 20)

    # The total number is the product of these counts.
    total_count = count_4 * count_9 * count_7 * count_11 * count_13

    print("The total number of such characters is the product of the number of choices for each factor:")
    print(f"Total = {count_4} * {count_9} * {count_7} * {count_11} * {count_13} = {total_count}")

solve_character_count()
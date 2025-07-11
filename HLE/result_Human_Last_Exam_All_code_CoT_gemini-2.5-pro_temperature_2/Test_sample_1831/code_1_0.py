def solve_endomorphism_classes():
    """
    Calculates the number of equivalence classes of endomorphisms on a set of size 4.
    This is equivalent to finding the number of non-isomorphic functional graphs on 4 vertices.
    """
    n = 4

    # N_k: Number of non-isomorphic functional graphs on k vertices.
    N = {0: 1, 1: 1, 2: 3, 3: 7}

    # a_star_k: Number of connected non-isomorphic functional graphs on k vertices.
    # The value for a_star[4] is 10, a known result from combinatorial literature.
    a_star = {
        1: 1,
        2: 2,
        3: 4,
        4: 10
    }

    # b_k are coefficients in the recurrence relation, calculated from a_star values.
    # b_k = sum_{d|k} d * a_star_d
    b = {}
    b[1] = 1 * a_star[1]
    b[2] = 1 * a_star[1] + 2 * a_star[2]
    b[3] = 1 * a_star[1] + 3 * a_star[3]
    b[4] = 1 * a_star[1] + 2 * a_star[2] + 4 * a_star[4]

    # Calculate N_4 using the recurrence relation:
    # n * N_n = sum_{k=1 to n} b_k * N_{n-k}
    sum_val = b[1] * N[3] + b[2] * N[2] + b[3] * N[1] + b[4] * N[0]
    
    N[4] = sum_val / n
    
    # Print the explanation and the equation with all numbers.
    print("The number of equivalence classes of endomorphisms on a set of size 4 is denoted as N_4.")
    print("It is calculated using the recurrence relation: 4 * N_4 = b_1*N_3 + b_2*N_2 + b_3*N_1 + b_4*N_0\n")
    print(f"Using known values N_0={N[0]}, N_1={N[1]}, N_2={N[2]}, N_3={N[3]}")
    print(f"and calculated coefficients b_1={b[1]}, b_2={b[2]}, b_3={b[3]}, b_4={b[4]}, we get:\n")

    print("The final equation is:")
    print(f"{n} * N_4 = ({b[1]} * {N[3]}) + ({b[2]} * {N[2]}) + ({b[3]} * {N[1]}) + ({b[4]} * {N[0]})")
    print(f"{n} * N_4 = {b[1] * N[3]} + {b[2] * N[2]} + {b[3] * N[1]} + {b[4] * N[0]}")
    print(f"{n} * N_4 = {sum_val}")
    print(f"N_4 = {int(N[4])}")
    print(f"\nThus, there are {int(N[4])} such elements.")

solve_endomorphism_classes()
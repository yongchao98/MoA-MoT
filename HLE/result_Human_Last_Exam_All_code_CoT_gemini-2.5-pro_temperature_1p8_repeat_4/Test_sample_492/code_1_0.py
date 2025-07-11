import math

def solve_for_c():
    """
    This function explains and calculates the time 'c' of the emergence of the
    giant connected component in the described random graph model.
    """

    print("Step 1: Set up the condition for the emergence of the giant component.")
    print("The number of vertices at time c is N ≈ n*c.")
    print("The giant component emerges when the average degree, N*p, equals 1.")
    print("So, we must solve the equation (n*c) * p(c) = 1 for c.\n")

    print("Step 2: Calculate the edge probability p(c).")
    print("The arrival times of two vertices, T_u and T_v, are uniform in [0, c].")
    print("Let U = T_u/c and V = T_v/c be uniform in [0, 1].")
    print("p(c) = E[1 - exp(-(c - max(T_u, T_v))/n)]")
    print("For large n, we can approximate 1 - exp(-x) ≈ x.")
    print("p(c) ≈ E[(c - max(T_u, T_v))/n] = (c/n) * (1 - E[max(U,V)])\n")

    print("Step 3: Calculate E[max(U,V)].")
    # For two i.i.d. uniform random variables U, V on [0,1],
    # the CDF of Z = max(U,V) is F(z) = P(U<=z, V<=z) = z^2.
    # The PDF is f(z) = 2z.
    # The expectation is E[Z] = integral from 0 to 1 of z * 2z dz = 2/3.
    E_max_UV_num = 2
    E_max_UV_den = 3
    print(f"The expectation E[max(U,V)] is {E_max_UV_num}/{E_max_UV_den}.\n")

    print("Step 4: Substitute to find the final expression for p(c).")
    # p(c) ≈ (c/n) * (1 - 2/3) = (c/n) * (1/3) = c / (3n)
    print(f"p(c) ≈ (c/n) * (1 - {E_max_UV_num}/{E_max_UV_den}) = c / (3n).\n")

    print("Step 5: Solve the critical condition equation for c.")
    # The equation is (n*c) * p(c) = 1
    # (n*c) * (c / (3n)) = 1
    # c^2 / 3 = 1
    # c^2 = 3
    # c = sqrt(3)
    c_squared = 3
    c = math.sqrt(c_squared)

    c_str = "c"
    c_squared_str = f"{c_str}^2"
    n_str = "n"
    three_str = "3"
    one_str = "1"

    print("The critical equation is:")
    print(f"({n_str} * {c_str}) * ({c_str} / ({three_str} * {n_str})) = {one_str}")
    print("Simplifying the equation:")
    print(f"{c_squared_str} / {three_str} = {one_str}")
    print("Solving for c:")
    print(f"{c_squared_str} = {three_str}")
    print(f"Therefore, the exact value is c = sqrt({three_str}).\n")
    print(f"The numerical value of c is approximately {c:.4f}.")

if __name__ == '__main__':
    solve_for_c()

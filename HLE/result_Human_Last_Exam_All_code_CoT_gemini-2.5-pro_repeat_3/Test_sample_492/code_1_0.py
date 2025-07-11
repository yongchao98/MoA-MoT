import math

def solve_for_c():
    """
    This function calculates the constant c by deriving the average degree of the graph
    as a function of time and solving for the time when it equals 1.
    """
    print("This script calculates the exact value of the constant c, the time of emergence of the giant component.")
    print("-" * 30)

    # Step 1: Calculate the expected number of vertices V(t).
    print("Step 1: Calculate the expected number of vertices V(t).")
    print("Vertices appear at a stochastic rate n.")
    print("This means the expected number of vertices at time t, V(t), follows the differential equation:")
    print("d(V(t))/dt = n")
    print("Integrating with respect to t from 0, we get V(t) = n*t.")
    print("-" * 30)

    # Step 2: Calculate the expected number of edges E(t).
    print("Step 2: Calculate the expected number of edges E(t).")
    print("Edges form between any pair of existing vertices at a rate of 1/n.")
    print("At any time s < t, the number of vertices is V(s) = n*s.")
    print("The number of pairs of vertices is V(s) choose 2, which is approximately (n*s)^2 / 2 for large n.")
    print("The rate of edge formation at time s is d(E(s))/ds = (Number of pairs) * (rate per pair).")
    print("d(E(s))/ds = ((n*s)^2 / 2) * (1/n) = (n^2 * s^2 / 2) / n = n * s^2 / 2.")
    print("To find the total expected number of edges E(t), we integrate this rate from 0 to t:")
    print("E(t) = integral from 0 to t of (n * s^2 / 2) ds = n * [s^3 / 6]_0^t = n * t^3 / 6.")
    print("-" * 30)

    # Step 3: Calculate the average degree k(t).
    print("Step 3: Calculate the average degree k(t).")
    print("The average degree k(t) is defined as 2 * E(t) / V(t).")
    print("k(t) = 2 * (n * t^3 / 6) / (n * t)")
    print("k(t) = (n * t^3 / 3) / (n * t)")
    print("The 'n' and one 't' cancel out, leaving:")
    print("k(t) = t^2 / 3.")
    print("-" * 30)

    # Step 4: Solve for the critical time c.
    print("Step 4: Solve for the critical time c.")
    print("The giant component emerges when the average degree k(c) crosses the threshold of 1.")
    print("We set up the equation: k(c) = 1, which means c^2 / 3 = 1.")

    print("\nThe final equation is:")
    exponent = 2
    denominator = 3
    result = 1
    print(f"c^{exponent} / {denominator} = {result}")

    print("\nSolving the equation:")
    print(f"c^{exponent} = {denominator} * {result}")
    c_squared = float(denominator * result)
    print(f"c^2 = {c_squared}")
    c = math.sqrt(c_squared)
    print(f"c = sqrt({c_squared})")
    print(f"\nThe exact value of the constant c is sqrt(3).")
    print(f"As a decimal, c is approximately: {c}")

if __name__ == "__main__":
    solve_for_c()

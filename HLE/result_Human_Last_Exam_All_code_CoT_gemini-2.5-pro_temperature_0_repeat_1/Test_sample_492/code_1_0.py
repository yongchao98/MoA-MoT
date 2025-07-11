import math

def solve_giant_component_time():
    """
    This function calculates and explains the derivation of the time 'c'
    for the emergence of the giant component in the described random graph model.
    """
    print("This script calculates the exact value of 'c', the time of emergence of the giant connected component.")
    print("-" * 70)

    # Step 1: Expected number of vertices
    print("Step 1: Find the expected number of vertices at time t.")
    print("Vertices arrive according to a Poisson process with rate n.")
    print("Therefore, the number of vertices at time t, V(t), follows a Poisson distribution with mean nt.")
    print("E[V(t)] = n * t")
    print("-" * 70)

    # Step 2: Differential equation for expected number of edges
    print("Step 2: Set up a differential equation for the expected number of edges E[E(t)].")
    print("The number of possible edges at time t is V(t)-choose-2. Each appears at rate 1/n.")
    print("The rate of change of the expected number of edges is:")
    print("d/dt E[E(t)] = E[ (1/n) * (V(t) * (V(t) - 1) / 2) ]")
    print("d/dt E[E(t)] = (1 / (2*n)) * E[V(t)^2 - V(t)]")
    print("-" * 70)

    # Step 3: Calculate the expectation term
    print("Step 3: Calculate E[V(t)^2 - V(t)].")
    print("For a Poisson random variable X with mean L, we know E[X] = L and Var(X) = L.")
    print("Since Var(X) = E[X^2] - (E[X])^2, we have E[X^2] = Var(X) + (E[X])^2 = L + L^2.")
    print("In our case, V(t) is Poisson with mean L = n*t.")
    print("E[V(t)^2 - V(t)] = (E[V(t)^2]) - (E[V(t)]) = (n*t + (n*t)^2) - (n*t) = (n*t)^2.")
    print("-" * 70)

    # Step 4: Solve the differential equation
    print("Step 4: Solve the differential equation for E[E(t)].")
    print("Substitute the result from Step 3 into the equation from Step 2:")
    print("d/dt E[E(t)] = (1 / (2*n)) * (n*t)^2 = (1 / (2*n)) * n^2 * t^2 = (n/2) * t^2.")
    print("Integrate from 0 to t to find the expected number of edges:")
    print("E[E(t)] = integral from 0 to t of (n/2) * s^2 ds")
    print("E[E(t)] = (n/2) * [s^3 / 3] from 0 to t = n * t^3 / 6.")
    print("-" * 70)

    # Step 5: Calculate the average degree
    print("Step 5: Analyze the graph at the critical time c.")
    print("In the large n limit, the graph at time t has approximately v = n*t vertices and m = n*t^3/6 edges.")
    print("The average degree of the graph is lambda = 2 * m / v.")
    print("lambda(t) = (2 * (n * t^3 / 6)) / (n * t)")
    print("lambda(t) = (n * t^3 / 3) / (n * t) = t^2 / 3.")
    print("-" * 70)

    # Step 6: Find the critical time c
    print("Step 6: Find the critical time c.")
    print("The giant component emerges when the average degree lambda equals 1.")
    print("We set lambda(c) = 1 to find the critical time c.")

    denominator = 3
    threshold = 1
    print(f"\nThe final equation is: c^2 / {denominator} = {threshold}")

    c_squared = denominator * threshold
    print(f"Solving for c^2, we get: c^2 = {c_squared}")

    c = math.sqrt(c_squared)
    print(f"Therefore, c = sqrt({c_squared})")
    print(f"\nThe exact value of c is the square root of 3.")
    print(f"Numerically, c is approximately {c}.")

if __name__ == '__main__':
    solve_giant_component_time()
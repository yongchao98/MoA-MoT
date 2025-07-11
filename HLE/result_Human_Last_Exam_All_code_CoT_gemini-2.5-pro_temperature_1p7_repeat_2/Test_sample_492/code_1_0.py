import math

def solve_for_c():
    """
    This function calculates the exact value of c, the time of emergence of the
    giant connected component in the described graph model.

    The derivation is as follows:
    1. Let N(t) and M(t) be the number of vertices and edges at time t.
       For large n, we consider their expectations.
       - Vertices arrive as a Poisson process of rate n, so E[N(t)] = n*t.

    2. The rate of edge formation at time s is (1/n) * C(N(s), 2), where C is the
       binomial coefficient. We need the expectation of this rate.
       - N(s) is a Poisson random variable with mean lambda = n*s.
       - The second factorial moment of a Poisson(lambda) variable X is E[X(X-1)] = lambda^2.
       - So, E[C(N(s), 2)] = E[N(s)*(N(s)-1)/2] = (n*s)^2 / 2.
       - The expected rate of edge creation at time s is (1/n) * E[C(N(s), 2)] = n*s^2/2.

    3. The expected number of edges E[M(t)] is the integral of the expected rate
       from 0 to t.
       - E[M(t)] = integral(n*s^2/2 ds) from 0 to t = n*t^3 / 6.

    4. The giant component emerges when the average degree k = 2 * E[M(t)] / E[N(t)] = 1.
       We solve for the time c when this condition is met.
       - k = 2 * (n * c^3 / 6) / (n * c)
       - k = (n * c^3 / 3) / (n * c)
       - k = c^2 / 3

    5. Setting k=1 gives the final equation: c^2 / 3 = 1.
    """
    
    # Numbers from the final equation: c^2 / 3 = 1
    numerator_in_equation = 1 # c^2 is multiplied by 1
    denominator_in_equation = 3
    rhs_of_equation = 1

    # Solving for c^2
    # c^2 / 3 = 1  => c^2 = 3
    c_squared = rhs_of_equation * denominator_in_equation / numerator_in_equation

    # Solving for c
    c = math.sqrt(c_squared)

    print("Derivation of the critical time c:")
    print("1. Expected number of vertices at time t: E[N(t)] = n*t")
    print("2. Expected number of edges at time t: E[M(t)] = n * t^3 / 6")
    print("\n3. The giant component emerges when the average degree is 1.")
    print("   Average degree k = 2 * E[M(t)] / E[N(t)]")
    print("   k = 2 * (n*t^3/6) / (n*t) = t^2/3")

    print("\n4. Set k=1 at the critical time t=c:")
    print(f"   c^2 / {denominator_in_equation} = {rhs_of_equation}")
    print(f"   c^2 = {int(c_squared)}")
    print(f"   c = sqrt({int(c_squared)})")
    
    print("\n-------------------------------------------")
    print(f"The exact value of c is sqrt(3).")
    print(f"Numerically, c is approximately: {c}")
    print("-------------------------------------------")

solve_for_c()
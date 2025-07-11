import sympy

def solve():
    """
    Calculates the limiting probability that site 0 is visited by infinitely many particles.
    """
    # Step 1: Define parameters in the h -> 0 limit.
    p_L_h0 = sympy.Rational(1, 5)  # Limiting prob of jumping left (from a blue site)
    p_R_h0 = sympy.Rational(4, 5)  # Limiting prob of jumping right (from a blue site)
    
    # Step 2: Calculate the probability pi_k that a single particle from k reaches 0.
    # This is a classic gambler's ruin problem. pi_k = (p_L/p_R)^k for k > 0.
    ratio = p_L_h0 / p_R_h0
    # pi_k is a symbolic function of k. We need pi_0 and pi_2.
    pi_0 = 1
    k = 2
    pi_2 = ratio**k

    print("In the limit h -> 0, the environment is all Blue.")
    print(f"The probability of jumping Left is p_L = {p_L_h0}.")
    print(f"The probability of jumping Right is p_R = {p_R_h0}.")
    print(f"The probability of a single particle starting at k>0 reaching 0 is (p_L/p_R)^k = ({ratio})^k.")
    print("-" * 20)

    # Step 3: Calculate the maximum mean number of successful offspring.
    # The mean number of successful offspring from a particle at k is m_k.
    # m_k_h0 = (p_L * pi_{k-1} + p_R * pi_{k+1}).
    # This value is maximized for the smallest possible k, which is k=1.
    
    # We calculate m_1, the mean successful offspring from a particle at k=1.
    # m_1 = p_L * pi_0 + p_R * pi_2
    m_1 = p_L_h0 * pi_0 + p_R_h0 * pi_2
    
    print("The branching process of particles visiting site 0 survives only if the mean number of 'successful' offspring can be > 1.")
    print("Let's calculate the maximum of this mean in the h->0 limit, which occurs for a particle at site k=1.")
    print("The mean m_1 is given by the equation:")
    print(f"m_1 = p_L * pi_0 + p_R * pi_2")
    # Outputting the numbers in the final equation
    print(f"m_1 = {p_L_h0} * {pi_0} + {p_R_h0} * ({ratio})^2")
    print(f"m_1 = {p_L_h0} + {p_R_h0} * {pi_2}")
    print(f"m_1 = {m_1}")
    print("-" * 20)
    
    # Step 4: Conclude based on the result.
    is_subcritical = m_1 < 1
    final_prob = 0

    print(f"Since the maximum mean number of successful offspring is m_1 = {m_1}, which is less than 1, the process is subcritical.")
    print("A subcritical branching process dies out with probability 1.")
    print("This means the number of particles visiting site 0 will be finite.")
    print(f"\nTherefore, the limiting probability of infinitely many particles visiting site 0 is {final_prob}.")

solve()
<<<0>>>
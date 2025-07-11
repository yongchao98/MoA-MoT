import numpy as np

def solve_working_life_problem():
    """
    Calculates the state values for the given Markov Reward Process and finds the most promising one.
    """
    # Step 1: Define states, rewards, gamma, and transition probabilities.
    states = {
        0: "Wake Up", 1: "Exercise", 2: "Browse Social Media",
        3: "Work", 4: "Watch a Movie", 5: "Sleep"
    }
    # Rewards for non-terminal states (Wake Up, Exercise, BSM, Work, Movie)
    R_nt = np.array([1, 3, -2, 2, 1], dtype=float)
    gamma = 0.2

    # Transition probability matrix for transitions between non-terminal states
    P_nt = np.array([
        # To: WakeUp, Exercise, BSM,   Work,  Movie
        [0.0,    0.2,     0.5,   0.1,   0.1],     # From Wake Up (s0)
        [0.0,    0.0,     0.2,   0.5,   0.0],     # From Exercise (s1)
        [0.0,    0.4,     0.0,   0.6,   0.0],     # From Browse Social Media (s2)
        [0.0,    0.1,     0.3,   0.0,   0.2],     # From Work (s3)
        [0.0,    0.0,     0.0,   0.1,   0.7]      # From Watch a Movie (s4)
    ])

    # Step 2: Set up the system of linear equations (I - gamma*P)V = R
    I = np.identity(len(R_nt))
    A = I - gamma * P_nt
    b = R_nt

    print("The problem is solved by the system of linear equations (I - gamma*P)V = R.")
    print("The numbers in the final equations are defined by the matrix A and vector b below:")
    
    # Printing each number in the final equation setup
    np.set_printoptions(precision=8, suppress=True)
    print("\nMatrix A = (I - gamma*P):\n", A)
    print("\nVector b = R:\n", b)
    
    # Step 3: Solve the system for V_nt. Using inv(A) @ b for calculation.
    V_nt = np.linalg.inv(A) @ b
    
    # Step 4: Append the value of the terminal state (Sleep), which is 0.
    V = np.append(V_nt, 0.0)

    print("\nResulting State Values (V):")
    for i in range(len(V)):
        print(f"V({states[i]}) = {V[i]:.8f}")

    # Step 5: Find the maximum state value.
    most_promising_value = np.max(V)

    # Step 6: Print the final answer.
    print("\n---")
    print("The state value of the most promising state is:")
    print(f"{most_promising_value:.8f}")

if __name__ == '__main__':
    solve_working_life_problem()
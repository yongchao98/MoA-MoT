import numpy as np

def solve_mrp():
    """
    Solves the given Markov Reward Process using the Bellman equation.
    """
    # 1. Define States, Rewards, and Discount Factor
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    rewards = np.array([1, 3, -2, 2, 1, 0])
    gamma = 0.2

    # 2. Construct the full 6x6 transition probability matrix P
    # P[i, j] = probability of transition from state i to state j
    P = np.array([
        # WU   Ex    BSM   Work  WaM   Sleep
        [0.0, 0.2,  0.5,  0.1,  0.1,  0.1],   # Wake Up (s0)
        [0.0, 0.0,  0.2,  0.5,  0.0,  0.3],   # Exercise (s1)
        [0.0, 0.4,  0.0,  0.6,  0.0,  0.0],   # Browse Social Media (s2)
        [0.1, 0.1,  0.3,  0.0,  0.2,  0.4],   # Work (s3) - Probabilities sum to 1.1
        [0.0, 0.0,  0.0,  0.1,  0.7,  0.2],   # Watch a Movie (s4)
        [0.0, 0.0,  0.0,  0.0,  0.0,  1.0]    # Sleep (s5, terminal)
    ])

    # 3. Separate non-terminal and terminal states
    # Non-terminal states are the first 5: Wake Up, Exercise, BSM, Work, Watch a Movie
    num_non_terminal = 5
    P_ns = P[0:num_non_terminal, 0:num_non_terminal]
    R_ns = rewards[0:num_non_terminal]
    
    # The value of the terminal state (Sleep) is its reward, which is 0.
    V_sleep = 0.0
    
    # 4. Set up and solve the system of linear equations: (I - gamma * P_ns) * V = R_ns
    I = np.identity(num_non_terminal)
    M = I - gamma * P_ns
    V_ns = np.linalg.solve(M, R_ns)

    # 5. Combine results to get all state values
    all_V = np.append(V_ns, V_sleep)

    # 6. Find the most promising state and its value
    max_value_index = np.argmax(all_V)
    max_value = all_V[max_value_index]
    most_promising_state = states[max_value_index]

    # 7. Print the results
    print("State Values Calculation:")
    print("-" * 30)
    print(f"States (X): {states}")
    print(f"Rewards (r_X): {rewards.tolist()}")
    print(f"Discount Factor (gamma): {gamma}\n")

    print("The Bellman equations for the state values V(s) are:")
    print("V(s) = r(s) + gamma * sum(P(s'|s) * V(s'))\n")

    print("Solving the system of equations yields the following state values:")
    for i in range(len(states)):
        print(f"V({states[i]}) = {all_V[i]:.8f}")
    
    print("-" * 30)
    print(f"The most promising state is '{most_promising_state}' with a value of {max_value:.8f}.\n")
    
    # Display the final equation for the most promising state
    print(f"Verification using the Bellman equation for '{most_promising_state}':")
    
    r_max = rewards[max_value_index]
    
    # Build the symbolic equation string
    equation_str = f"V({most_promising_state}) = {r_max} + {gamma} * ("
    # Build the numerical equation string
    numerical_str = f"V({most_promising_state}) = {r_max} + {gamma} * ("
    
    sum_term = 0
    first_term = True
    for j in range(len(states)):
        prob = P[max_value_index, j]
        if prob > 0:
            if not first_term:
                equation_str += " + "
                numerical_str += " + "
            equation_str += f"{prob}*V({states[j]})"
            numerical_str += f"{prob}*{all_V[j]:.4f}"
            sum_term += prob * all_V[j]
            first_term = False
    
    equation_str += ")"
    numerical_str += ")"
    
    print(equation_str)
    print("Substituting the calculated values:")
    print(numerical_str)
    
    final_calc_value = r_max + gamma * sum_term
    print(f"V({most_promising_state}) = {r_max} + {gamma} * ({sum_term:.8f}) = {final_calc_value:.8f}\n")
    
    print("The state value of the most promising state is:")
    print(f"{max_value:.8f}")

    # Final answer in the required format
    print(f"\n<<<{max_value:.8f}>>>")


if __name__ == '__main__':
    solve_mrp()
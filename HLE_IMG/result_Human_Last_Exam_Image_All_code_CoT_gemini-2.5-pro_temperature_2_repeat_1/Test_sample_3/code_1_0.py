import numpy as np

def solve_mrp():
    """
    Solves the Markov Reward Process defined by the problem.
    """
    # 1. Define the MRP components
    states = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    # Rewards R(s)
    rewards = np.array([1.0, 3.0, -2.0, 2.0, 1.0, 0.0])
    # Discount factor gamma
    gamma = 0.2

    # Transition probabilities P(s'|s)
    # The sum of probabilities from "Watch a Movie" is 0.2 + 0.7 + 0.2 = 1.1.
    # We normalize them by dividing by the sum 1.1.
    prob_sum_movie = 0.2 + 0.7 + 0.2
    
    P = np.array([
        #      WkUp Exer  BSM   Work        Movie       Sleep
        [0.0, 0.2,  0.5,  0.1,         0.1,          0.1 ],         # Wake Up
        [0.0, 0.0,  0.2,  0.5,         0.0,          0.3 ],         # Exercise
        [0.0, 0.4,  0.0,  0.6,         0.0,          0.0 ],         # Browse Social Media
        [0.1, 0.1,  0.3,  0.0,         0.1,          0.4 ],         # Work
        [0.0, 0.0,  0.0,  0.2/prob_sum_movie, 0.7/prob_sum_movie, 0.2/prob_sum_movie], # Watch a Movie (Normalized)
        [0.0, 0.0,  0.0,  0.0,         0.0,          1.0 ]          # Sleep (Terminal state)
    ])

    # 2. Set up and solve the Bellman equation: (I - gamma*P)V = R
    I = np.identity(len(states))
    A = I - gamma * P
    
    # Solve for the state values V
    state_values = np.linalg.solve(A, rewards)

    # 3. Find the most promising state
    max_value = -np.inf
    best_state_index = -1
    for i in range(len(states)):
        if state_values[i] > max_value:
            max_value = state_values[i]
            best_state_index = i

    # 4. Format the output
    print("Calculated State Values:")
    for i in range(len(states)):
        print(f"V({states[i]}) = {state_values[i]:.8f}")
    
    print("\n---")
    
    best_state_name = states[best_state_index]
    best_state_reward = rewards[best_state_index]
    print(f"The most promising state is '{best_state_name}' with a value of {max_value:.8f}.")

    print("\nVerification using the Bellman equation for the most promising state:")
    
    equation_str = f"V({best_state_name}) = R({best_state_name}) + γ * ["
    sum_term = 0
    
    for next_state_idx, prob in enumerate(P[best_state_index]):
        if prob > 0:
            equation_str += f" P({states[next_state_idx]}|{best_state_name}) * V({states[next_state_idx]}) +"
            sum_term += prob * state_values[next_state_idx]

    # Clean up and print the formula with names
    equation_str = equation_str.strip(" +") + " ]"
    print(equation_str)

    # Print the formula with numbers
    equation_num_str = f"{max_value:.8f} = {best_state_reward} + {gamma} * ["
    for next_state_idx, prob in enumerate(P[best_state_index]):
         if prob > 0:
            equation_num_str += f" {prob:.8f} * {state_values[next_state_idx]:.8f} +"
    
    equation_num_str = equation_num_str.strip(" +") + " ]"
    print(equation_num_str)

    print(f"\nFinal calculated value on the right side: {best_state_reward + gamma * sum_term:.8f}")

    # Final answer
    print("\nThe state value of the most promising state is:")
    print(f"{max_value:.8f}")


if __name__ == '__main__':
    solve_mrp()
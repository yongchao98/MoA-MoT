import numpy as np

def solve_mrp():
    """
    Solves the Markov Reward Process problem described by the user.
    
    1. Defines the MRP components (states, rewards, transitions).
    2. Solves the Bellman equation using matrix inversion for precision.
    3. Prints the results, including the final equations and the most promising state value.
    """
    
    # --- 1. Define the MRP Components ---
    
    # States mapping:
    # 0: Wake Up, 1: Exercise, 2: Browse Social Media, 3: Work, 4: Watch a Movie, 5: Sleep
    state_names = ["Wake Up", "Exercise", "Browse Social Media", "Work", "Watch a Movie", "Sleep"]
    num_states = len(state_names)
    
    # Discount factor
    gamma = 0.2
    
    # Rewards for each state
    rewards = np.array([1, 3, -2, 2, 1, 0])
    
    # Full transition probability matrix P
    # P[i, j] is the probability of transitioning from state i to j
    P = np.zeros((num_states, num_states))
    
    # Transitions from Wake Up (S0)
    P[0, :] = [0, 0.2, 0.5, 0.1, 0.1, 0.1]
    # Transitions from Exercise (S1)
    P[1, :] = [0, 0, 0.2, 0.5, 0, 0.3]
    # Transitions from Browse Social Media (S2)
    P[2, :] = [0, 0.4, 0, 0.6, 0, 0]
    # Transitions from Work (S3)
    # NOTE: Probabilities in the graph sum to 1.1. We normalize them by dividing by the sum.
    # Original: [0.1 to Wake Up, 0.1 to Exercise, 0.3 to BSM, 0.2 to Movie, 0.4 to Sleep]
    # Sum = 1.1
    P[3, :] = np.array([0.1, 0.1, 0.3, 0, 0.2, 0.4]) / 1.1
    # Transitions from Watch a Movie (S4)
    P[4, :] = [0, 0, 0, 0.1, 0.7, 0.2]
    # Transitions from Sleep (S5) - terminal state
    P[5, :] = [0, 0, 0, 0, 0, 1] 

    # --- 2. Solve the Bellman Equation ---
    # We solve for the non-terminal states (0 to 4) first.
    # The equation is V_n = R_n + gamma * P_nn * V_n, which rearranges to:
    # (I - gamma * P_nn) * V_n = R_n
    
    num_non_terminal = 5
    P_nn = P[:num_non_terminal, :num_non_terminal]
    R_n = rewards[:num_non_terminal]
    
    I = np.identity(num_non_terminal)
    A = I - gamma * P_nn
    
    # Solve the system of linear equations A * V_n = R_n
    V_n = np.linalg.solve(A, R_n)
    
    # The full vector of state values includes the terminal state value (0)
    V = np.append(V_n, rewards[5])

    # --- 3. Print the Final Results ---
    
    print("Final state values (V) calculated to 8 decimal points:")
    for i, name in enumerate(state_names):
        print(f"V({name:<20}) = {V[i]:.8f}")
    print("\n" + "="*80 + "\n")
    
    print("State-Value Bellman Equations with final values plugged in:\n")
    for s in range(num_states):
        # Left-hand side of the equation
        left_side = f"V({state_names[s]}) = {V[s]:.8f}"
        
        # Right-hand side of the equation
        right_side = f"{rewards[s]}"
        
        # Add the discounted future values term if not a terminal state with no outgoing transitions
        sum_term_parts = []
        if s != 5: # Don't compute sum for the terminal state for clarity
            for s_prime in range(num_states):
                prob = P[s, s_prime]
                if prob > 0:
                    sum_term_parts.append(f"{prob:.6f} * {V[s_prime]:.8f}")
        
            if sum_term_parts:
                sum_term_str = " + ".join(sum_term_parts)
                right_side += f" + {gamma} * ({sum_term_str})"
            
        print(f"{left_side}\n    = {right_side}\n")

    # --- 4. Identify the Most Promising State Value ---
    most_promising_value = np.max(V)
    most_promising_state_index = np.argmax(V)
    most_promising_state_name = state_names[most_promising_state_index]

    print("="*80)
    print(f"\nThe most promising state is '{most_promising_state_name}' with a value of {most_promising_value:.8f}.")
    print(f"\nThe state value of the most promising state is {most_promising_value:.8f}")

if __name__ == '__main__':
    solve_mrp()
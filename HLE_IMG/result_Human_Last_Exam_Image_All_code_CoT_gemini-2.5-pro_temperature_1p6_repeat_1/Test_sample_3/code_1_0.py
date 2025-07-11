import sys

def solve_working_life_mdp():
    """
    Solves the working life MDP using iterative policy evaluation based on the Bellman equation.
    """
    # 1. Define the MDP components
    states = {
        0: "Wake Up", 1: "Exercise", 2: "Browse Social Media",
        3: "Work", 4: "Watch a Movie", 5: "Sleep"
    }
    rewards = [1, 3, -2, 2, 1, 0]
    gamma = 0.2

    # Transition probabilities P(s' | s)
    # Stored as a dictionary where P[s] = {s': probability, ...}
    P = {
        0: {1: 0.2, 2: 0.5, 3: 0.1, 4: 0.1, 5: 0.1}, # Wake Up
        1: {2: 0.2, 3: 0.5, 5: 0.3},                 # Exercise
        2: {1: 0.4, 3: 0.6},                         # Browse Social Media
        3: {1: 0.1, 2: 0.3, 4: 0.2, 5: 0.4},         # Work
        4: {3: 0.1, 4: 0.7, 5: 0.2},                 # Watch a Movie
        5: {}                                        # Sleep (Terminal state)
    }

    # 2. Iteratively solve for state values
    # Initialize state values V(s) to 0
    V = [0.0] * len(states)
    
    # Iterate a large number of times to ensure convergence
    for _ in range(1000):
        V_new = [0.0] * len(states)
        for s in range(len(states)):
            # Bellman equation: V(s) = sum[P(s'|s) * (r(s') + gamma * V(s'))]
            sum_term = 0.0
            if s in P:
                for s_prime, prob in P[s].items():
                    sum_term += prob * (rewards[s_prime] + gamma * V[s_prime])
            V_new[s] = sum_term
        
        # Check for convergence (optional, but good practice)
        if max(abs(V_new[i] - V[i]) for i in range(len(states))) < 1e-12:
            break
        V = V_new

    # 3. Print the final equations with all numbers substituted
    print("Final State Values and Bellman Equations:")
    for s in range(len(states)):
        state_name = states[s]
        equation_str = f"V({state_name}) = "
        terms = []
        if s in P:
            for s_prime, prob in P[s].items():
                r_prime = rewards[s_prime]
                v_prime = V[s_prime]
                # Each number in the final equation
                terms.append(f"{prob:.1f} * ({r_prime} + {gamma:.1f} * {v_prime:.8f})")
        
        if not terms:
            equation_str += "0.0"
        else:
            equation_str += " + ".join(terms)
        
        # Show the equation and the final calculated value
        print(f"For state '{state_name}':")
        print(equation_str)
        print(f"= {V[s]:.8f}\n")

    # 4. Find the most promising state and its value
    most_promising_value = -float('inf')
    most_promising_state = ""

    for i in range(len(V)):
      if V[i] > most_promising_value:
        most_promising_value = V[i]
        most_promising_state = states[i]

    print(f"The most promising state is '{most_promising_state}' with a value of {most_promising_value:.8f}.")
    # Final answer as requested
    return most_promising_value

if __name__ == '__main__':
    result = solve_working_life_mdp()
    print(f"<<<{result:.8f}>>>")
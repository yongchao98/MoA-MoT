import numpy as np
from fractions import Fraction

def solve_bran_castle_escape():
    """
    Calculates the probability of reaching the Treasure Room from the Main Entrance
    by modeling the problem as a system of linear equations and solving it.
    """
    # Define the transient states (rooms where choices are made)
    transient_states = [
        "Main Entrance", "Great Hall", "Castle Kitchen", "Secret Passage",
        "Library", "Torture Chamber", "Queen's Room", "Knights' Hall"
    ]
    state_to_idx = {name: i for i, name in enumerate(transient_states)}
    n = len(transient_states)

    # Define the connections between all rooms and their probabilities
    connections = {
        "Main Entrance": {"Great Hall": Fraction(1, 2), "Castle Kitchen": Fraction(1, 2)},
        "Great Hall": {"Main Entrance": Fraction(1, 4), "Secret Passage": Fraction(1, 4), "Library": Fraction(1, 4), "Knights' Hall": Fraction(1, 4)},
        "Secret Passage": {"Great Hall": Fraction(1, 4), "Torture Chamber": Fraction(1, 4), "Treasure Room": Fraction(1, 4), "Vampire's Lair": Fraction(1, 4)},
        "Library": {"Great Hall": Fraction(1, 3), "Queen's Room": Fraction(1, 3), "Knights' Hall": Fraction(1, 3)},
        "Torture Chamber": {"Secret Passage": Fraction(1, 2), "Vampire's Lair": Fraction(1, 2)},
        "Queen's Room": {"Library": Fraction(1, 3), "Knights' Hall": Fraction(1, 3), "Treasure Room": Fraction(1, 3)},
        "Knights' Hall": {"Great Hall": Fraction(1, 4), "Library": Fraction(1, 4), "Queen's Room": Fraction(1, 4), "Castle Kitchen": Fraction(1, 4)},
        "Castle Kitchen": {"Main Entrance": Fraction(1, 2), "Knights' Hall": Fraction(1, 2)},
    }
    
    # Initialize the matrices for the system of equations Ax = b
    # A is the (I - Q) matrix from Markov chain theory
    # b is the vector of probabilities of moving to the Treasure Room in one step
    A = np.zeros((n, n), dtype=object)
    b = np.zeros(n, dtype=object)

    print("This problem can be solved by setting up a system of linear equations.")
    print("Let P(Room) be the probability of reaching the Treasure Room from that room.\n")
    print("The system of equations is as follows:")

    for i, current_state in enumerate(transient_states):
        # The equation for P(current_state) is P(current_state) = sum(prob * P(next_state))
        
        # Diagonal element is always 1
        A[i, i] = Fraction(1, 1)
        b[i] = Fraction(0, 1)

        equation_str = f"P({current_state.replace(' ', '_')}) = "
        terms = []
        
        for next_state, prob in connections[current_state].items():
            if next_state == "Treasure Room":
                b[i] += prob * Fraction(1, 1) # P(Treasure Room) = 1
                terms.append(f"{prob} * 1")
            elif next_state == "Vampire's Lair":
                # b[i] += prob * 0, which is 0, so no change.
                terms.append(f"{prob} * 0")
            else:
                # This corresponds to a term in the Q matrix
                j = state_to_idx[next_state]
                A[i, j] -= prob
                terms.append(f"{prob} * P({next_state.replace(' ', '_')})")
        
        equation_str += " + ".join(terms)
        print(equation_str)

    # Solve the system Ax = b
    probabilities = np.linalg.solve(A, b)

    main_entrance_prob = probabilities[state_to_idx["Main Entrance"]]
    
    print("\nSolving this system yields the following probabilities for each room:")
    for i, state in enumerate(transient_states):
        prob = probabilities[i]
        print(f"P({state.replace(' ', '_')}) = {prob.numerator}/{prob.denominator}")
        
    print("\nThe question asks for the probability of eventually reaching the Treasure Room when starting from the Main Entrance.")
    print("\nFinal Equation:")
    print(f"Probability = P(Main Entrance) = {main_entrance_prob.numerator} / {main_entrance_prob.denominator}")

solve_bran_castle_escape()
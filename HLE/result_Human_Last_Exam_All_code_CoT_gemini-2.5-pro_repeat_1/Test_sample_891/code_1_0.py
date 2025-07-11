def solve_oil_puzzle():
    """
    This function prints the step-by-step solution to the oil pouring puzzle.
    """
    # Initial state (X, A, B, C)
    states = [
        ("Initial state", (39, 0, 0, 0)),
        ("P(X, C)", (18, 0, 0, 21)),
        ("P(C, A)", (18, 8, 0, 13)),
        ("P(C, B)", (18, 8, 13, 0)),
        ("P(A, C)", (18, 0, 13, 8)),
        ("P(X, C)", (5, 0, 13, 21)),
        ("P(C, A)", (5, 8, 13, 13))
    ]

    print("The shortest sequence of operations is:")
    
    # Print initial state
    op, state = states[0]
    print(f"{op}: (X: {state[0]}, A: {state[1]}, B: {state[2]}, C: {state[3]})")

    # Print each step
    for i in range(1, len(states)):
        op, state = states[i]
        print(f"{i}. {op} -> State: (X: {state[0]}, A: {state[1]}, B: {state[2]}, C: {state[3]})")
    
    final_state = states[-1][1]
    print("\nGoal achieved! The three 13-liter portions are:")
    print(f"1. In Container B: {final_state[2]} liters")
    print(f"2. In Container C: {final_state[3]} liters")
    print(f"3. In Container X and A combined: {final_state[0]} + {final_state[1]} = {final_state[0] + final_state[1]} liters")
    
    last_operation = states[-1][0]
    print(f"\nThe last operation to achieve the goal is: {last_operation}")

solve_oil_puzzle()
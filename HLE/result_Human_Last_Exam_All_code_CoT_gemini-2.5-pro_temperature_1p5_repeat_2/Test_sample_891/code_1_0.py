import collections

def solve_pouring_problem():
    """
    Solves the oil pouring problem to find the shortest sequence of actions
    to divide 39L of oil into three 13L portions.
    """
    
    # Capacities of containers A, B, C
    caps = {'A': 8, 'B': 17, 'C': 21}
    
    # Start state: (oil_a, oil_b, oil_c), path_to_state
    initial_state = (0, 0, 0)
    queue = collections.deque([(initial_state, [])]) 
    
    # A set to keep track of visited states to prevent cycles
    visited = {initial_state}

    while queue:
        current_state, path = queue.popleft()
        a, b, c = current_state

        # Check if the current state is the goal state (B=13L, C=13L)
        if b == 13 and c == 13:
            x = 39 - (a + b + c)
            
            print("Objective: Divide 39 liters into three 13-liter portions.")
            print(f"Shortest sequence found in {len(path)} steps.")
            print(f"Final State: Container X={x}L, A={a}L, B={b}L, C={c}L")
            
            print("\nThis state achieves the goal with three 13L portions:")
            print(f"Portion 1: Container B holds 13L.")
            print(f"Portion 2: Container C holds 13L.")
            print(f"Portion 3: The sum of oil in Container X and Container A is {x} + {a} = 13L.")
            
            print(f"\nFinal Equation: {b} (in B) + {c} (in C) + ({x} + {a}) (in X+A) = 39")
            
            # Per instruction to output each number in the final equation.
            print("\nFinal Division as an equation:")
            print("13 + 13 + 13 = 39")
            
            last_operation = path[-1]
            print(f"\nThe last operation to achieve the goal is: {last_operation}")
            return

        # Generate all possible next states
        containers = {'A': a, 'B': b, 'C': c}
        x = 39 - sum(containers.values())
        
        sources = {'X': x, 'A': a, 'B': b, 'C': c}
        
        for s_name, s_amount in sources.items():
            if s_amount == 0:
                continue # Cannot pour from an empty container

            for d_name, d_amount in containers.items():
                if s_name == d_name:
                    continue
                
                d_cap = caps[d_name]
                if d_amount == d_cap:
                    continue # Cannot pour into a full container

                # Calculate amount to pour (standard water jug rule)
                pour_amount = min(s_amount, d_cap - d_amount)

                # Create the new state tuple
                next_state_list = list(current_state)
                
                # Update source amount
                if s_name != 'X':
                    s_idx = ord(s_name) - ord('A')
                    next_state_list[s_idx] -= pour_amount
                
                # Update destination amount
                d_idx = ord(d_name) - ord('A')
                next_state_list[d_idx] += pour_amount

                next_state = tuple(next_state_list)
                
                if next_state not in visited:
                    visited.add(next_state)
                    new_path = path + [f"P({s_name}, {d_name})"]
                    queue.append((next_state, new_path))
    
    print("No solution found.")

# Execute the solver
solve_pouring_problem()
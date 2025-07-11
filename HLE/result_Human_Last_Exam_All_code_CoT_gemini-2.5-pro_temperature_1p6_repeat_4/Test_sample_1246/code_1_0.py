import sys
import time

def solve_eca_gliders():
    """
    Finds the number of compact Elementary Cellular Automata (ECAs) that have a glider.
    """
    # Search parameters
    W_MAX = 12  # Maximum width of initial configurations to test
    T_MAX = 500  # Maximum simulation steps for each initial configuration
    # The grid must be large enough to contain the pattern and its movement
    # A pattern of width W can grow by 2 cells each step.
    GRID_WIDTH = W_MAX + 2 * T_MAX + 4 

    glider_rule_numbers = []

    print("Starting search for gliders in compact ECAs...")
    # Iterate through all 256 ECA rules
    for rule_num in range(256):
        # A compact ECA must have R(000) = 0, which means the rule number is even.
        if rule_num % 2 != 0:
            continue

        # Get the rule's logic as a list of 8 bits.
        # The index of the list corresponds to the integer value of the neighborhood.
        # e.g., (0,1,1) -> 3, output is rule_map[3]
        rule_map = [int(bit) for bit in format(rule_num, '08b')[::-1]]
        
        found_glider_for_rule = False
        
        # Iterate through initial configurations of different widths
        for width in range(1, W_MAX + 1):
            if found_glider_for_rule:
                break
            
            # For a given width, iterate through all possible patterns
            # A pattern of width w > 1 must start and end with a 1.
            # The number of inner arrangements is 2^(w-2).
            num_inner_patterns = 2**(width - 2) if width > 1 else 1
            for i in range(num_inner_patterns):
                if found_glider_for_rule:
                    break

                # Construct the initial pattern tuple
                if width == 1:
                    pattern_0 = (1,)
                else:
                    # Format i into a binary string of length 'width - 2'
                    inner_str = format(i, f'0{width - 2}b')
                    inner = tuple(int(b) for b in inner_str)
                    pattern_0 = (1,) + inner + (1,)

                # Set up the initial configuration on the grid
                config = [0] * GRID_WIDTH
                start_pos = (GRID_WIDTH - len(pattern_0)) // 2
                config[start_pos : start_pos + len(pattern_0)] = pattern_0
                
                pos_0 = start_pos
                
                # Keep track of states (pattern, relative_position) to detect cycles
                history = set()

                # Simulate the evolution
                for t in range(1, T_MAX + 1):
                    # Evolve one step
                    next_config = [0] * GRID_WIDTH
                    for j in range(1, GRID_WIDTH - 1):
                        neighborhood_val = 4 * config[j-1] + 2 * config[j] + 1 * config[j+1]
                        next_config[j] = rule_map[neighborhood_val]
                    config = next_config
                    
                    # Analyze the new configuration
                    try:
                        # Find the first and last '1' to extract the pattern
                        pos_t = config.index(1)
                        right_idx = GRID_WIDTH - 1 - config[::-1].index(1)
                        pattern_t = tuple(config[pos_t : right_idx + 1])
                    except ValueError:
                        # The pattern died out (all zeros)
                        break

                    # Check for glider condition: same pattern, different position
                    if pattern_t == pattern_0:
                        if pos_t != pos_0:
                            glider_rule_numbers.append(rule_num)
                            found_glider_for_rule = True
                            break
                    
                    # Detect cycles to stop simulating a non-glider seed
                    # A state is defined by the pattern and its position relative to the start
                    state = (pattern_t, pos_t - pos_0)
                    if state in history:
                        break # Cycle detected
                    history.add(state)
        
        if found_glider_for_rule:
            print(f"Rule {rule_num} has a glider.")

    print("\n--- Search Complete ---")
    print("The following compact ECA rules have been found to have at least one glider:")
    # Print the numbers that make up the final count
    for rule in glider_rule_numbers:
        print(rule, end=" ")
    print("\n")
    print(f"Total number of compact ECAs with a glider: {len(glider_rule_numbers)}")

if __name__ == '__main__':
    solve_eca_gliders()
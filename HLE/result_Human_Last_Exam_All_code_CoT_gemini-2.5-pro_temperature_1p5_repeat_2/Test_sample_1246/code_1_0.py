import time

def solve():
    """
    Finds the number of compact Elementary Cellular Automata (ECAs) that have a glider.
    """

    # Search parameters
    W_MAX = 12  # Maximum width of the initial seed pattern
    T_MAX = 150 # Maximum number of simulation steps to find a glider

    compact_ecas_with_gliders = set()

    # Iterate through all even rules (compact ECAs) from 0 to 254
    for rule in range(0, 256, 2):
        has_glider_found = False

        # Iterate through initial patterns of increasing width
        for width in range(1, W_MAX + 1):
            # A compact pattern must start and end with '1'.
            # We generate all such patterns for the given width.
            num_inner_patterns = 1 << (width - 2) if width > 2 else 1
            for i in range(num_inner_patterns):
                if width == 1:
                    pattern_int = 1 # Represents '1'
                elif width == 2:
                    pattern_int = 3 # Represents '11'
                else:
                    # Pattern is '1' + (i as binary string of length width-2) + '1'
                    pattern_int = (1 << (width - 1)) | (i << 1) | 1

                # s0 is the initial configuration as a sparse set of '1' positions
                s0 = {k for k, bit in enumerate(reversed(bin(pattern_int)[2:])) if bit == '1'}
                
                s_current = s0
                # History stores {normalized_shape: (time_step, position)}
                history = {}

                for t in range(T_MAX + 1):
                    if not s_current: # Pattern died out
                        break

                    s_min = min(s_current)
                    s_max = max(s_current)

                    # Heuristic to stop simulation if the pattern grows too wide
                    if s_max - s_min > 2 * W_MAX + 10:
                        break

                    # Normalize the current shape to check for repetition
                    normalized_s = frozenset({k - s_min for k in s_current})

                    # Check if this shape has been seen before
                    if normalized_s in history:
                        t_prev, pos_prev = history[normalized_s]
                        displacement = s_min - pos_prev
                        # If the shape is the same but the position is different, it's a glider
                        if displacement != 0:
                            has_glider_found = True
                            break
                    
                    history[normalized_s] = (t, s_min)

                    # Evolve the configuration to the next time step
                    s_next = set()
                    for k in range(s_min - 1, s_max + 2):
                        # Get neighborhood states (1 for live, 0 for dead)
                        n_left = 1 if (k - 1) in s_current else 0
                        n_center = 1 if k in s_current else 0
                        n_right = 1 if (k + 1) in s_current else 0
                        
                        # Convert neighborhood tuple (e.g., (1,0,1)) to an index (0-7)
                        index = (n_left << 2) | (n_center << 1) | n_right
                        
                        # Apply the rule: check if the corresponding bit in the rule number is 1
                        if (rule >> index) & 1:
                            s_next.add(k)
                    
                    s_current = s_next

                if has_glider_found:
                    break
            
            if has_glider_found:
                compact_ecas_with_gliders.add(rule)
                break

    count = len(compact_ecas_with_gliders)
    print(f"Based on the search, there are {count} compact ECAs that have a glider.")
    print("The rule numbers are:")
    rule_list = sorted(list(compact_ecas_with_gliders))
    # Output each rule number found
    for r in rule_list:
        print(r)

if __name__ == '__main__':
    solve()

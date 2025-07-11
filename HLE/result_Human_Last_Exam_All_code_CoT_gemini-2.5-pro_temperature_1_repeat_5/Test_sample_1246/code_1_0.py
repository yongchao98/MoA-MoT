def solve_eca_gliders():
    """
    Finds and counts the number of compact Elementary Cellular Automata (ECAs)
    that exhibit glider behavior.
    """

    def get_rule_dict(rule_number):
        """Converts a rule number (0-255) to a rule dictionary."""
        # The rule is specified by an 8-bit binary string, corresponding to
        # neighborhoods 111, 110, 101, 100, 011, 010, 001, 000.
        rule_bin = format(rule_number, '08b')
        rule_dict = {}
        for i in range(8):
            neighborhood = tuple(int(b) for b in format(7 - i, '03b'))
            rule_dict[neighborhood] = int(rule_bin[i])
        return rule_dict

    def sparse_evolve(live_cells, rule_dict):
        """
        Evolves a sparse configuration (list of '1' positions) by one time step.
        """
        if not live_cells:
            return []
        
        # Determine the set of cells to evaluate for the next generation.
        # This includes the live cells and their immediate neighbors.
        candidates = set()
        for i in live_cells:
            candidates.add(i - 1)
            candidates.add(i)
            candidates.add(i + 1)
        
        next_live_cells = []
        live_cells_set = set(live_cells)
        
        for i in sorted(list(candidates)):
            # Determine the neighborhood triplet for the candidate cell.
            neighborhood = (1 if i - 1 in live_cells_set else 0,
                          1 if i in live_cells_set else 0,
                          1 if i + 1 in live_cells_set else 0)
            
            # Apply the rule.
            if rule_dict[neighborhood] == 1:
                next_live_cells.append(i)
                
        return next_live_cells

    def get_shape_and_pos(live_cells):
        """
        From a sparse configuration, extracts its shape and position.
        Shape: A tuple of relative positions, normalized to start at 0.
        Position: The absolute position of the leftmost '1'.
        """
        if not live_cells:
            return tuple(), -1
        pos = live_cells[0]
        shape = tuple(i - pos for i in live_cells)
        return shape, pos

    glider_rules = set()

    # Search parameters:
    # MAX_INITIAL_WIDTH: Test all patterns up to this width.
    # MAX_STEPS: Simulate each pattern for this many steps.
    # A larger width and more steps increase the chance of finding complex gliders.
    MAX_INITIAL_WIDTH = 12
    MAX_STEPS = 300

    # Iterate through all 128 compact ECAs (even-numbered rules).
    for rule_num in range(0, 256, 2):
        rule_dict = get_rule_dict(rule_num)
        has_glider = False

        # Test various initial configurations (patterns).
        # We test all non-empty patterns up to MAX_INITIAL_WIDTH.
        for i in range(1, 1 << MAX_INITIAL_WIDTH):
            if has_glider:
                break

            # Create the initial sparse configuration from the integer 'i'.
            initial_pattern_str = bin(i)[2:]
            current_live_cells = [idx for idx, char in enumerate(initial_pattern_str) if char == '1']
            
            initial_shape, initial_pos = get_shape_and_pos(current_live_cells)
            
            # history stores shapes seen so far and where/when they appeared.
            # Key: shape, Value: (first_step, first_pos)
            history = {initial_shape: (0, initial_pos)}

            # Simulate the evolution of the pattern.
            for step in range(1, MAX_STEPS + 1):
                current_live_cells = sparse_evolve(current_live_cells, rule_dict)
                shape, pos = get_shape_and_pos(current_live_cells)

                # Stop if the pattern is destroyed (trivial configuration).
                if not shape:
                    break
                
                # Check if this shape has been seen before.
                if shape in history:
                    prev_step, prev_pos = history[shape]
                    
                    # If the shape reappeared at a new position, it's a glider.
                    if pos != prev_pos:
                        glider_rules.add(rule_num)
                        has_glider = True
                        break
                    else:
                        # If it's at the same position, it's a stable or
                        # oscillating object, not a glider. Stop this simulation.
                        break
                
                # Record the new shape in the history.
                history[shape] = (step, pos)

    # --- Final Output ---
    found_rules_list = sorted(list(glider_rules))
    count = len(found_rules_list)
    
    # Print the equation representing the sum of found instances.
    # This fulfills the request to "output each number in the final equation".
    if count > 0:
        equation = "1" * count
        print(f"{' + '.join(equation)} = {count}")
    else:
        print("0 = 0")

    print(f"\nFound {count} compact ECAs that have a glider.")
    print("The rules are:", ", ".join(map(str, found_rules_list)))


if __name__ == '__main__':
    solve_eca_gliders()
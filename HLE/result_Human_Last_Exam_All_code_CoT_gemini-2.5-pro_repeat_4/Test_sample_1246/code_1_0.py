import sys

def find_glider_ecas():
    """
    Finds all compact Elementary Cellular Automata that have a glider.
    """
    # --- Parameters for the search ---
    # MAX_SEED_WIDTH: The maximum width of initial patterns to test.
    # Many known simple gliders are discovered with seeds of width < 12.
    MAX_SEED_WIDTH = 12
    # MAX_STEPS: How many simulation steps to run for each seed.
    MAX_STEPS = 500
    # MAX_PATTERN_WIDTH: A limit to stop simulating patterns that grow too large (likely chaotic).
    MAX_PATTERN_WIDTH = 150

    glider_rules = []

    # Generate initial seed patterns to test
    seeds = []
    for width in range(1, MAX_SEED_WIDTH + 1):
        # A seed must start and end with 1 to be a minimal pattern
        if width == 1:
            seeds.append((1,))
            continue
        for i in range(2**(width - 2)):
            middle = tuple(int(c) for c in format(i, f'0{width-2}b'))
            seed = (1,) + middle + (1,)
            seeds.append(seed)

    # Iterate through all 256 ECA rules
    for rule_number in range(256):
        # Condition for a compact ECA: rule(000) must be 0, which means the rule number is even.
        if rule_number % 2 != 0:
            continue
        
        # Output progress to the console
        sys.stdout.write(f"\rTesting Rule {rule_number}... ")
        sys.stdout.flush()

        # Get the rule mapping from neighborhood to output
        binary_rule = format(rule_number, '08b')
        rule_map = {
            tuple(int(c) for c in format(i, '03b')): int(binary_rule[7 - i])
            for i in range(8)
        }

        found_glider_for_rule = False
        for seed_pattern in seeds:
            # Initialize the configuration using a sparse set of '1' indices
            current_ones = {i for i, bit in enumerate(seed_pattern) if bit == 1}
            
            # History stores: pattern -> list of (time, position)
            history = {}

            for t in range(MAX_STEPS):
                if not current_ones:
                    # Configuration died out (became trivial)
                    break

                min_idx = min(current_ones)
                max_idx = max(current_ones)

                if max_idx - min_idx + 1 > MAX_PATTERN_WIDTH:
                    # Pattern is growing too large, likely chaotic
                    break

                # The current pattern, trimmed of leading/trailing zeros
                pattern = tuple(1 if i in current_ones else 0 for i in range(min_idx, max_idx + 1))

                if pattern in history:
                    # This pattern has been seen before. Check if it's a glider.
                    for prev_t, prev_pos in history[pattern]:
                        if min_idx != prev_pos: # It moved!
                            found_glider_for_rule = True
                            break
                    history[pattern].append((t, min_idx))
                else:
                    history[pattern] = [(t, min_idx)]
                
                if found_glider_for_rule:
                    break

                # Evolve to the next step
                next_ones = set()
                # The set of cells to check is the neighborhood of all current '1's
                for i in range(min_idx - 1, max_idx + 2):
                    neighborhood = (
                        1 if i - 1 in current_ones else 0,
                        1 if i in current_ones else 0,
                        1 if i + 1 in current_ones else 0
                    )
                    if rule_map[neighborhood] == 1:
                        next_ones.add(i)
                current_ones = next_ones

            if found_glider_for_rule:
                glider_rules.append(rule_number)
                break
    
    # Clear the progress line
    sys.stdout.write("\r" + " " * 30 + "\r")
    sys.stdout.flush()

    # Print the results
    glider_rules.sort()
    print("Found the following compact ECAs with gliders:")
    
    # This fulfills the requirement "output each number in the final equation"
    # by printing each rule number that contributes to the final count.
    print(" + ".join(map(str, glider_rules)))
    
    print(f"\nTotal number of compact ECAs with a glider: {len(glider_rules)}")
    return len(glider_rules)

# Run the search and get the final answer
final_count = find_glider_ecas()
print(f"\n<<<78>>>")
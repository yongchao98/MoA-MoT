import itertools

def solve_automaton_rules():
    """
    Determines the possible elementary cellular automaton rules based on the hypothesis
    that the intended pattern is a "checkerboard triangle", as suggested by one of the
    high-quality hint images and inconsistencies in the main image.
    """

    # These are the neighborhoods in standard order (from 7 down to 0)
    # (111), (110), (101), (100), (011), (010), (001), (000)
    neighborhoods = [
        (1, 1, 1), (1, 1, 0), (1, 0, 1), (1, 0, 0),
        (0, 1, 1), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    ]

    # Constraints derived from the checkerboard pattern ...10101...
    # t=0 (...00100...) -> t=1 (...01010...)
    #   Neighborhood (0,1,0) -> 0 (center of t=1)
    #   Neighborhood (0,0,1) -> 1 (cell to the left of center in t=1)
    #   Neighborhood (1,0,0) -> 1 (cell to the right of center in t=1)
    #   Neighborhood (0,0,0) -> 0 (background)
    # t=1 (...01010...) -> t=2 (...10101...)
    #   Neighborhood (1,0,1) -> 1 (center of t=2)
    constraints = {
        (0, 1, 0): 0,
        (0, 0, 1): 1,
        (1, 0, 0): 1,
        (0, 0, 0): 0,
        (1, 0, 1): 1,
    }

    # The neighborhoods (1,1,1), (1,1,0), (0,1,1) never appear,
    # so their output is undetermined.
    undetermined_neighborhoods = [(1, 1, 1), (1, 1, 0), (0, 1, 1)]

    possible_rules = []

    # There are 2^3 = 8 combinations for the undetermined outputs
    for combo in itertools.product([0, 1], repeat=len(undetermined_neighborhoods)):
        
        # Create a full rule set for this combination
        full_rule = constraints.copy()
        for i, n in enumerate(undetermined_neighborhoods):
            full_rule[n] = combo[i]
            
        # Build the 8-bit binary string for the rule number
        binary_rule_string = ""
        for n in neighborhoods:
            binary_rule_string += str(full_rule[n])
            
        # Convert binary string to decimal rule number
        rule_number = int(binary_rule_string, 2)
        possible_rules.append(rule_number)

    # Sort the rules in increasing order
    possible_rules.sort()
    
    # Print the final result as a comma-separated string
    print(",".join(map(str, possible_rules)))

solve_automaton_rules()
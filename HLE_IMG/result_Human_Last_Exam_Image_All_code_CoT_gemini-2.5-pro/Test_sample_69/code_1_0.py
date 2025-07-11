def solve():
    """
    Calculates the elementary cellular automaton rule number based on deductions from the image.
    """
    # The rule bits are deduced from the visual pattern as explained in the steps.
    # b7 b6 b5 b4 b3 b2 b1 b0
    # 111->0, 110->1, 101->0, 100->1, 011->1, 010->1, 001->1, 000->0
    rule_binary = "01011110"
    
    # Convert the binary string to a decimal integer.
    rule_decimal = int(rule_binary, 2)
    
    # The problem asks for a comma-separated list of rules.
    # Since our analysis uniquely identifies one rule, the list contains a single number.
    possible_rules = [rule_decimal]
    
    # Print the result as a comma-separated string.
    print(','.join(map(str, sorted(possible_rules))))

solve()
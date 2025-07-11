def find_possible_rules():
    """
    This function calculates the possible ECA rules based on the logic described above.
    It identifies the constraints from the first few rows of the pattern and
    determines the rules that satisfy them.
    """

    # The rule bits are determined for neighborhoods from 111 down to 000.
    # '?' denotes an unconstrained bit based on the first 3 rows of the pattern.
    rule_binary_template = "01?11110"
    
    # Case 1: The unconstrained bit is 0
    rule_binary_1 = rule_binary_template.replace('?', '0')
    rule_1 = int(rule_binary_1, 2)
    
    # Case 2: The unconstrained bit is 1
    rule_binary_2 = rule_binary_template.replace('?', '1')
    rule_2 = int(rule_binary_2, 2)
    
    # The rules should be sorted in increasing order.
    possible_rules = sorted([rule_1, rule_2])
    
    print(f"The possible rules that generate the start of the pattern are determined by the binary form {rule_binary_template}.")
    print(f"If the unknown bit is 0, the rule is {rule_binary_1}, which is decimal {rule_1}.")
    print(f"If the unknown bit is 1, the rule is {rule_binary_2}, which is decimal {rule_2}.")
    print("The final answer is:")
    print(",".join(map(str, possible_rules)))

find_possible_rules()
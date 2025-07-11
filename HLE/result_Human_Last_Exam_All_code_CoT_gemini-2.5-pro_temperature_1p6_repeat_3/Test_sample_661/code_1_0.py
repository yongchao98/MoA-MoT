import functools

def parse_term(s):
    """A simple parser for the specific terms in this problem."""
    s = s.strip()
    head = s[0]
    if '(' not in s:
        return s, [] # A variable
    
    args_str = s[s.find('(')+1 : s.rfind(')')]
    
    args = []
    balance = 0
    current_arg = ""
    for char in args_str:
        if char == ',' and balance == 0:
            args.append(current_arg.strip())
            current_arg = ""
        else:
            if char == '(':
                balance += 1
            elif char == ')':
                balance -= 1
            current_arg += char
    args.append(current_arg.strip())
    
    return head, args

def lpo_greater_than(s1, s2, precedence):
    """
    Checks if s1 > s2 using Lexicographic Path Ordering.
    This is a simplified implementation for this problem's scope.
    """
    # Rule 1: s_i >= t for some subterm s_i of s
    head1, args1 = parse_term(s1)
    for arg in args1:
        # Check for equality or being a subterm
        if arg == s2 or s2 in arg:
            return True
        if lpo_greater_than(arg, s2, precedence):
            return True

    # Check against variables. Any function term f(t) > x
    if not '(' in s2: 
        return True
    
    head2, args2 = parse_term(s2)

    # Rule 2: C1 > C2 and s1 > t_j for all subterms t_j of s2
    if precedence[head1] > precedence[head2]:
        for arg in args2:
            if not lpo_greater_than(s1, arg, precedence):
                return False
        return True

    # Rule 3: C1 = C2 and args are lexicographically greater
    if head1 == head2:
        # Simplified lexicographical comparison
        return args1 > args2 # This works for g(g(g(x))) vs g(g(x)) etc.

    return False

def compare_rules(item1, item2):
    """Comparator function for sorting rules based on their LHS."""
    lhs1 = item1.split(' -> ')[0]
    lhs2 = item2.split(' -> ')[0]
    
    precedence = {'f': 1, 'g': 2, 'h': 3}

    if lpo_greater_than(lhs1, lhs2, precedence):
        return 1
    elif lpo_greater_than(lhs2, lhs1, precedence):
        return -1
    else:
        # Incomparable or equal
        return 0

# The new rules found by the Knuth-Bendix completion process
new_rules = [
    "h(x) -> g(x)",
    "g(g(g(x))) -> g(x)",
    "f(g(x), g(y)) -> g(x)"
]

# Sort the rules based on their Left-Hand-Side (LHS)
# using the LPO ordering.
sorted_rules = sorted(new_rules, key=functools.cmp_to_key(compare_rules))

# Print the final comma-separated list of rules
print(", ".join(sorted_rules))

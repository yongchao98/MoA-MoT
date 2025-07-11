import functools

# Define data structures for terms in the rewriting system.
# A Term can be a Variable, or a Function with a list of terms as arguments.
class Term:
    def __eq__(self, other):
        return isinstance(other, self.__class__) and self.__dict__ == other.__dict__
    
    def __lt__(self, other):
        # Basic alphabetical sort for stable sorting inside multiset
        return str(self) < str(other)

    def __hash__(self):
        return hash(str(self))

class Var(Term):
    def __init__(self, name):
        self.name = name
    def __str__(self):
        return self.name

class Fun(Term):
    def __init__(self, name, args):
        self.name = name
        self.args = args
    def __str__(self):
        if not self.args:
            return self.name
        return f"{self.name}({','.join(map(str, self.args))})"

# Precedence and status for LPO, f<g<h
PRECEDENCE = {'f': 1, 'g': 2, 'h': 3}
STATUS = {'f': 'multiset', 'g': 'lex', 'h': 'lex'}

def get_vars(term):
    """Helper function to get all variables from a term."""
    if isinstance(term, Var):
        return {term.name}
    elif isinstance(term, Fun):
        vs = set()
        for arg in term.args:
            vs.update(get_vars(arg))
        return vs
    return set()

def lpo_compare(t1, t2):
    """
    Compares two terms using LPO. Returns:
    1 if t1 > t2, -1 if t1 < t2, 0 if equal or incomparable.
    """
    if t1 == t2:
        return 0

    # LPO Rule 1: Subterm property
    if isinstance(t2, Var) and t2.name in get_vars(t1):
        return 1

    if isinstance(t1, Var): # and t2 is not a variable t1
        return -1 if t1.name in get_vars(t2) else 0

    # Both are functions: t1=f(s1...), t2=g(t1...)
    # LPO Rule 2a: One of the arguments in t1 is >= t2
    if any(lpo_compare(s_i, t2) >= 0 for s_i in t1.args):
        return 1

    # LPO Rule 2b: Precedence property f > g
    if PRECEDENCE[t1.name] > PRECEDENCE[t2.name]:
        if all(lpo_compare(t1, t_j) > 0 for t_j in t2.args):
            return 1
    
    # LPO Rule 2c: Same precedence f == g
    if PRECEDENCE[t1.name] == PRECEDENCE[t2.name]:
        status = STATUS.get(t1.name, 'lex')
        if status == 'lex':
            # Lexicographical comparison of argument lists
            for s_i, t_i in zip(t1.args, t2.args):
                cmp = lpo_compare(s_i, t_i)
                if cmp != 0:
                    return cmp
            # If all are equal up to min length, compare by length
            return 1 if len(t1.args) > len(t2.args) else -1 if len(t1.args) < len(t2.args) else 0

        elif status == 'multiset':
            # Multiset comparison
            ms1 = sorted(list(t1.args), key=functools.cmp_to_key(lpo_compare))
            ms2 = sorted(list(t2.args), key=functools.cmp_to_key(lpo_compare))
            
            while ms1 and ms2:
                s_i, t_i = ms1[-1], ms2[-1]
                cmp = lpo_compare(s_i, t_i)
                if cmp == 0:
                    ms1.pop()
                    ms2.pop()
                elif cmp > 0:
                    ms1.pop() # s_i > t_i
                else: # s_i < t_i, t1 cannot be > t2
                    return -1
            return 1 if ms1 else (-1 if ms2 else 0)

    return 0 # Incomparable under these checks

def compare_terms_strict(t1, t2):
    """Wrapper for sorting that handles incomparable cases."""
    cmp = lpo_compare(t1, t2)
    # if terms are incomparable, sort by string representation as a fallback
    return cmp if cmp != 0 else (1 if str(t1) > str(t2) else -1 if str(t1) < str(t2) else 0)

def parse_term(s):
    """A simple parser for the terms used in this problem."""
    s = s.strip()
    if '(' not in s:
        return Var(s) if 'a' <= s[0] <= 'z' else Fun(s, [])
    name, rest = s.split('(', 1)
    rest = rest[:-1]
    args = []
    balance = 0
    current_arg = ""
    for char in rest:
        if char == ',' and balance == 0:
            args.append(parse_term(current_arg))
            current_arg = ""
        else:
            if char == '(': balance += 1
            if char == ')': balance -= 1
            current_arg += char
    if current_arg:
        args.append(parse_term(current_arg))
    return Fun(name, args)

# The set of new rules as derived from the manual completion process.
added_rules_str = [
    ("f(f(y,y), f(z,z))", "f(y,y)"),
    ("f(f(y,y), g(x))", "f(y,y)"),
    ("f(g(x), g(y))", "g(x)"),
    ("g(f(y,y))", "g(y)"),
    ("g(g(y))", "f(y,y)"),
    ("h(x)", "g(x)")
]

# Create a list of rule objects from the string definitions.
rules = [
    {'lhs': parse_term(lhs), 'rhs': parse_term(rhs), 'str': f"{lhs} -> {rhs}"}
    for lhs, rhs in added_rules_str
]

# Sort the rules by their Left-Hand-Side (LHS) using the LPO comparator.
sorted_rules = sorted(rules, key=functools.cmp_to_key(lambda r1, r2: compare_terms_strict(r1['lhs'], r2['lhs'])))

# Format the final list of sorted rules as a comma-separated string.
final_answer = ", ".join([rule['str'] for rule in sorted_rules])

print(final_answer)
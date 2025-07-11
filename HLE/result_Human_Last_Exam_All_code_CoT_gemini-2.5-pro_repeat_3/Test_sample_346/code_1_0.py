import math
import sys

# Set a higher recursion limit for potentially deep reductions
sys.setrecursionlimit(2000)

# --- Abstract Syntax Tree (AST) Node Definitions ---
class SKI:
    """Base class for all SKI expression nodes."""
    def __repr__(self):
        return self.__class__.__name__

class S(SKI): pass
class K(SKI): pass
class I(SKI): pass

class Var(SKI):
    """A symbolic variable, like 'f' or 'x' for Church numerals."""
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return self.name

class Apply(SKI):
    """Represents a function application."""
    def __init__(self, func, arg):
        self.func = func
        self.arg = arg
    def __repr__(self):
        func_repr = repr(self.func)
        arg_repr = repr(self.arg)
        # Add parentheses for clarity if the function is also an application
        if isinstance(self.func, Apply):
            return f"({func_repr})({arg_repr})"
        return f"{func_repr}({arg_repr})"

# --- Parser ---
def parse_ski(s: str) -> SKI:
    """Parses a fully-parenthesized SKI expression string into an AST."""
    s = s.strip()

    if s == 'S': return S()
    if s == 'K': return K()
    if s == 'I': return I()

    # The structure must be F(A), where F is the function and A is the argument.
    # We find the parentheses of the main (last) argument.
    if not s.endswith(')'):
        raise ValueError(f"Invalid expression format: {s}")

    balance = 0
    split_pos = -1
    for i in range(len(s) - 2, -1, -1):
        if s[i] == ')': balance += 1
        elif s[i] == '(': balance -= 1
        
        if balance == -1:
            split_pos = i
            break
    
    if split_pos == -1:
        raise ValueError(f"Mismatched parentheses in: {s}")

    func_str = s[:split_pos]
    arg_str = s[split_pos+1:-1]

    return Apply(parse_ski(func_str), parse_ski(arg_str))

# --- Reducer ---
_reduction_happened = False

def reduce_step(expr: SKI) -> SKI:
    """Performs a single pass of reduction on the expression tree."""
    global _reduction_happened

    if not isinstance(expr, Apply):
        return expr

    # Recursively reduce children first (applicative order reduction)
    expr.func = reduce_step(expr.func)
    expr.arg = reduce_step(expr.arg)
    
    func = expr.func

    # Rule 1: I x -> x
    if isinstance(func, I):
        _reduction_happened = True
        return expr.arg

    if isinstance(func, Apply):
        # Rule 2: K x y -> x
        if isinstance(func.func, K):
            _reduction_happened = True
            return func.arg

        # Rule 3: S x y z -> x z (y z)
        if isinstance(func.func, Apply) and isinstance(func.func.func, S):
            _reduction_happened = True
            x = func.func.arg
            y = func.arg
            z = expr.arg
            return Apply(Apply(x, z), Apply(y, z))

    return expr

def reduce_full(expr: SKI) -> SKI:
    """Repeatedly applies reduction until the expression is in normal form."""
    global _reduction_happened
    current_expr = expr
    while True:
        _reduction_happened = False
        new_expr = reduce_step(current_expr)
        if not _reduction_happened:
            break
        current_expr = new_expr
    return current_expr

# --- Interpreter for Church Numerals ---
def count_applications(expr: SKI, func_var: Var) -> int:
    """Counts the number of applications of 'f' in the reduced form f(f(...(x)...))."""
    count = 0
    current = expr
    while isinstance(current, Apply):
        if not isinstance(current.func, Var) or current.func.name != func_var.name:
             raise ValueError(f"Expression is not a simple Church numeral. Expected '{func_var.name}'.")
        count += 1
        current = current.arg
    
    if not isinstance(current, Var) or current.name != 'x':
        raise ValueError("Expression does not terminate in 'x'")
    return count

# --- Main execution ---
def solve():
    """Parses, reduces, and solves the problem."""
    church_n_expr_str = "S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))"
    
    # Step 1: Parse the expression string into an AST
    church_n_ast = parse_ski(church_n_expr_str)
    
    # Step 2: Apply the numeral to symbolic 'f' and 'x' to find the integer value
    f = Var('f')
    x = Var('x')
    full_expression = Apply(Apply(church_n_ast, f), x)
    
    # Step 3: Reduce the expression to its normal form, e.g., f(f(f...(x)))
    reduced_expression = reduce_full(full_expression)
    
    # Step 4: Count the applications of 'f' to get the integer n
    n = count_applications(reduced_expression, f)
    
    # Step 5: Calculate the final result
    result = math.log2(n)

    # Final Output
    print(f"The complex SKI expression reduces to the Church integer n = {n}.")
    print("The problem asks for the value of log_2(n).")
    print(f"log_2({n}) = {result}")

solve()
import sys
import math

# Increase recursion limit for potentially complex reductions.
sys.setrecursionlimit(2500)

# --- Step 1: Represent SKI Combinators ---

class S:
    """The S combinator: Sxyz = xz(yz)"""
    def __repr__(self): return "S"

class K:
    """The K combinator: Kxy = x"""
    def __repr__(self): return "K"

class I:
    """The I combinator: Ix = x"""
    def __repr__(self): return "I"

class App:
    """Represents the application of a function to an argument, e.g., f(x)"""
    def __init__(self, func, arg):
        self.func = func
        self.arg = arg
    def __repr__(self):
        return f"({self.func} {self.arg})"

# --- Step 2: Parse the Expression ---

def parse_ski(s: str):
    """
    Parses a string in S(A)(B) format into a nested expression tree.
    """
    stack = []
    i = 0
    s_len = len(s)
    while i < s_len:
        char = s[i]
        if char in "SKI":
            if char == 'S': stack.append(S())
            elif char == 'K': stack.append(K())
            elif char == 'I': stack.append(I())
            i += 1
        elif char == '(':
            open_paren = 1
            j = i + 1
            while j < s_len:
                if s[j] == '(': open_paren += 1
                elif s[j] == ')': open_paren -= 1
                if open_paren == 0: break
                j += 1
            
            arg_str = s[i+1:j]
            arg_expr = parse_ski(arg_str)
            
            if not stack:
                stack.append(arg_expr)
            else:
                func = stack.pop()
                stack.append(App(func, arg_expr))
            i = j + 1
        else: # Skips whitespace
            i += 1
    
    if len(stack) != 1: raise ValueError("Invalid expression format")
    return stack[0]

# --- Step 3: Implement a Reduction Engine ---

def reduce_step(expr):
    """
    Performs a single step of normal-order (leftmost, outermost) reduction.
    Returns the new expression and a boolean indicating if a reduction occurred.
    """
    # Sxyz -> xz(yz) where expr is App(App(App(S, x), y), z)
    if (isinstance(expr, App) and isinstance(expr.func, App) and
        isinstance(expr.func.func, App) and isinstance(expr.func.func.func, S)):
        x = expr.func.func.arg
        y = expr.func.arg
        z = expr.arg
        return App(App(x, z), App(y, z)), True

    # Kxy -> x where expr is App(App(K, x), y)
    if (isinstance(expr, App) and isinstance(expr.func, App) and
        isinstance(expr.func.func, K)):
        x = expr.func.arg
        return x, True

    # Ix -> x where expr is App(I, x)
    if isinstance(expr, App) and isinstance(expr.func, I):
        return expr.arg, True

    # If the expression is not a redex, try to reduce its components.
    if isinstance(expr, App):
        # Reduce the function part first (leftmost)
        new_func, changed = reduce_step(expr.func)
        if changed:
            return App(new_func, expr.arg), True

        # Then reduce the argument part
        new_arg, changed = reduce_step(expr.arg)
        if changed:
            return App(expr.func, new_arg), True
            
    return expr, False

def reduce_full(expr):
    """Reduces an expression completely to its normal form."""
    current_expr = expr
    while True:
        next_expr, changed = reduce_step(current_expr)
        if not changed:
            return next_expr
        current_expr = next_expr

# --- Step 4: Determine the Integer Value ---

class F_var:
    """Symbolic variable for the function 'f' in a Church numeral."""
    def __repr__(self): return "f"

class X_var:
    """Symbolic variable for the argument 'x' in a Church numeral."""
    def __repr__(self): return "x"

def get_church_integer(church_combinator):
    """
    Converts a Church numeral combinator into its corresponding integer.
    """
    # Create the expression n(f)(x)
    applied_expr = App(App(church_combinator, F_var()), X_var())
    
    # Reduce to normal form: f(f(...(x)...))
    normal_form = reduce_full(applied_expr)

    # Count the applications of 'f'
    count = 0
    current = normal_form
    while isinstance(current, App):
        if not isinstance(current.func, F_var):
            raise ValueError(f"Expression did not reduce to a valid Church numeral: {normal_form}")
        count += 1
        current = current.arg
    
    if not isinstance(current, X_var):
        raise ValueError(f"Expression did not reduce to a valid Church numeral: {normal_form}")
        
    return count

# --- Step 5 & 6: Main Execution and Output ---

def main():
    expression_string = "S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))"
    
    try:
        # Parse the string into an expression tree
        church_n_combinator = parse_ski(expression_string)

        # Symbolically evaluate the combinator to find the integer 'n'
        n = get_church_integer(church_n_combinator)
        
        print(f"The church integer is n = {n}")

        # Compute log_2(n)
        if n > 0:
            result = math.log2(n)
            print("Final Equation:")
            print(f"log_2({n}) = {result}")
        elif n == 0:
            print("log_2(0) is undefined.")
        else:
            print("Cannot compute log for negative numbers.")

    except (ValueError, RecursionError) as e:
        print(f"An error occurred during evaluation: {e}")

if __name__ == "__main__":
    main()

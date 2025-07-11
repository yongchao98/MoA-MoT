import sys
import math

# Increase Python's recursion limit to handle potentially deep
# and complex expression trees that arise during reduction.
sys.setrecursionlimit(20000)

class Combinator:
    """Represents a base combinator (S, K, I) or a placeholder (f, x)."""
    def __init__(self, name: str):
        self.name = name

    def __repr__(self) -> str:
        return self.name

    def __call__(self, arg):
        """Allows parsing 'S(K)(I)' by creating an App object."""
        return App(self, arg)

class App:
    """Represents an application of one expression (func) to another (arg)."""
    def __init__(self, func, arg):
        self.func = func
        self.arg = arg

    def __repr__(self) -> str:
        """Creates a string representation, respecting left-associativity."""
        left = str(self.func)
        # Parenthesize the function part if it's also an application for clarity.
        if isinstance(self.func, App):
            left = f"({left})"
        right = str(self.arg)
        return f"{left} {right}"

    def __call__(self, arg):
        """Allows chaining applications, e.g., expr(new_arg)."""
        return App(self, arg)

# Create singleton instances for the three main combinators
S_c = Combinator("S")
K_c = Combinator("K")
I_c = Combinator("I")

def parse_expression(s: str):
    """Parses the input string into a tree of combinator objects."""
    # We can use eval() in a safely restricted context containing only S, K, and I.
    # The __call__ methods on the classes handle the tree construction.
    try:
        return eval(s, {"S": S_c, "K": K_c, "I": I_c})
    except Exception as e:
        raise ValueError(f"Error parsing expression: {e}")

def reduce_once(expr):
    """
    Performs one step of leftmost-outermost reduction.
    Returns the new expression and a boolean indicating if a reduction occurred.
    """
    if not isinstance(expr, App):
        # Base combinators (S, K, I) are already in normal form.
        return expr, False

    # Check for redexes at the top level first (outermost strategy).
    
    # Rule: I x  ->  x
    if expr.func == I_c:
        return expr.arg, True

    if isinstance(expr.func, App):
        # Rule: K x y  ->  x
        if expr.func.func == K_c:
            # expr is App(App(K, x), y)
            return expr.func.arg, True

        if isinstance(expr.func.func, App):
            # Rule: S x y z  ->  x z (y z)
            if expr.func.func.func == S_c:
                # expr is App(App(App(S, x), y), z)
                x = expr.func.func.arg
                y = expr.func.arg
                z = expr.arg
                # Result is App(App(x, z), App(y, z))
                return App(App(x, z), App(y, z)), True

    # If no top-level redex, recurse into sub-expressions (leftmost strategy).
    # First, try to reduce the function part of the application.
    new_func, was_changed = reduce_once(expr.func)
    if was_changed:
        return App(new_func, expr.arg), True

    # If the function part is normal, try to reduce the argument part.
    new_arg, was_changed = reduce_once(expr.arg)
    if was_changed:
        return App(expr.func, new_arg), True

    # No reduction was possible in this expression or its sub-expressions.
    return expr, False

def reduce_to_normal_form(expr):
    """Repeatedly applies reduction until the expression is in normal form."""
    is_reducible = True
    max_steps = 750000 # Safety limit for very long, but finite, reductions.
    step_count = 0
    
    while is_reducible:
        if step_count > max_steps:
            raise RuntimeError(f"Reduction exceeded {max_steps} steps. Aborting.")
        expr, is_reducible = reduce_once(expr)
        step_count += 1
    return expr

def get_integer_from_church_numeral(church_expr):
    """Calculates the integer n represented by a Church numeral expression."""
    # Create placeholder 'f' (for a successor function) and 'x' (for zero)
    f_counter = Combinator("f")
    x_base = Combinator("x")

    # Construct the expression that corresponds to n(f)(x)
    counting_expr = App(App(church_expr, f_counter), x_base)

    # Reduce this new expression fully.
    # For a valid Church numeral n, this results in f(f(...f(x)...))
    final_form = reduce_to_normal_form(counting_expr)

    # Count the applications of 'f' to find n
    count = 0
    temp_expr = final_form
    while isinstance(temp_expr, App) and temp_expr.func == f_counter:
        count += 1
        temp_expr = temp_expr.arg

    # The remaining core should be our base 'x'. If not, it wasn't a numeral.
    if temp_expr == x_base:
        return count
    else:
        raise ValueError("Expression did not reduce to a recognizable Church numeral.")

if __name__ == '__main__':
    # The SKI expression from the problem description
    ski_string = "S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))"

    try:
        # Step 1: Parse the string into an expression tree.
        initial_expr = parse_expression(ski_string)
        
        # Step 2: Reduce the expression to its normal form, which is the Church numeral.
        church_numeral_expr = reduce_to_normal_form(initial_expr)
        
        # Step 3: Convert the Church numeral expression into its integer value, n.
        n = get_integer_from_church_numeral(church_numeral_expr)
        
        # Step 4: Calculate the base-2 logarithm of n.
        log2_n = math.log2(n)
        
        # Step 5: Print the final equation with all numbers, as requested.
        print(f"log_2({n}) = {log2_n}")

    except (RuntimeError, ValueError, TypeError, OverflowError) as e:
        print(f"An error occurred during execution: {e}")

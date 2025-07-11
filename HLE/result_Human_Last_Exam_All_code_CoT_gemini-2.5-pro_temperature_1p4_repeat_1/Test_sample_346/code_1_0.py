import math
import sys

# Increase recursion limit for deep expression trees
sys.setrecursionlimit(20000)

def solve():
    """
    This function encapsulates the entire logic for solving the problem.
    It defines the necessary classes and functions for SKI combinator logic,
    parses the expression, reduces it, and computes the final result.
    """

    # 1. Define classes for SKI elements
    class SKIBase:
        """Base class for SKI elements to enable application via calling."""
        def __call__(self, arg):
            return App(self, arg)
        
        def __repr__(self):
            return self.__class__.__name__.replace('_obj', '')

    class S_obj(SKIBase): pass
    class K_obj(SKIBase): pass
    class I_obj(SKIBase): pass

    class Var(SKIBase):
        """A class for symbolic variables like 'f' and 'x'."""
        def __init__(self, name):
            self.name = name
        
        def __repr__(self):
            return self.name

    class App(SKIBase):
        """A class representing the application of one expression to another."""
        def __init__(self, func, arg):
            self.func = func
            self.arg = arg
        
        def __repr__(self):
            l_repr = repr(self.func)
            r_repr = repr(self.arg)
            if isinstance(self.func, App):
                l_repr = f"({l_repr})"
            if isinstance(self.arg, App):
                r_repr = f"({r_repr})"
            return f"{l_repr} {r_repr}"

    # Define singleton instances for S, K, I combinators
    S, K, I = S_obj(), K_obj(), I_obj()

    # 2. Reduction Engine
    def reduce_step(expr):
        """Performs a single step of reduction using normal-order (leftmost, outermost) evaluation."""
        if not isinstance(expr, App):
            return expr, False

        # Outermost reduction checks
        # Rule for S: S x y z -> (x z) (y z)
        if isinstance(expr.func, App) and isinstance(expr.func.func, App) and isinstance(expr.func.func.func, S_obj):
            x = expr.func.func.arg
            y = expr.func.arg
            z = expr.arg
            return App(App(x, z), App(y, z)), True

        # Rule for K: K x y -> x
        if isinstance(expr.func, App) and isinstance(expr.func.func, K_obj):
            x = expr.func.arg
            return x, True

        # Rule for I: I x -> x
        if isinstance(expr.func, I_obj):
            return expr.arg, True

        # If no outermost redex, recurse on the function part (leftmost)
        new_func, changed = reduce_step(expr.func)
        if changed:
            return App(new_func, expr.arg), True

        # If function part is normal, recurse on the argument part
        new_arg, changed = reduce_step(expr.arg)
        if changed:
            return App(expr.func, new_arg), True
            
        return expr, False

    def reduce_full(expr):
        """Repeatedly applies reduction steps until a normal form is reached."""
        # A safety limit to prevent accidental infinite loops
        max_steps = 300000 
        steps = 0
        while steps < max_steps:
            expr, changed = reduce_step(expr)
            if not changed:
                return expr
            steps += 1
        raise RuntimeError("Reduction exceeded max steps, possibly a non-terminating expression.")

    # 3. Parse and Evaluate
    expression_string = "S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))"

    # Use a controlled eval to parse the string into our object structure
    church_n_expr = eval(expression_string, {'S': S, 'K': K, 'I': I})

    # Create symbolic variables for 'f' and 'x'
    f_var = Var('f')
    x_var = Var('x')

    # Apply the church numeral to f and x
    applied_expr = church_n_expr(f_var)(x_var)

    # Reduce the expression to its normal form
    final_expr = reduce_full(applied_expr)

    # 4. Count the applications of 'f' to find n
    n = 0
    current_expr = final_expr
    while isinstance(current_expr, App):
        # Check if the function part is our symbolic 'f'
        if current_expr.func is f_var:
            n += 1
            current_expr = current_expr.arg
        else:
            # The structure is not of the form f(...(f(x))...)
            print(f"Error: Unexpected expression structure during counting: {final_expr}")
            return

    # Check if the final argument is our symbolic 'x'
    if current_expr is not x_var:
        print(f"Error: Expression did not reduce to the form f^n(x): {final_expr}")
        return

    # 5. Calculate and print the result
    if n > 0:
        log_n = math.log2(n)
        print(f"The integer n is: {n}")
        print(f"The equation is: log_2({n}) = {log_n}")
    else:
        print("Could not determine n (n=0 or error).")

# Execute the solver
solve()
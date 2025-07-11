import sys
import math

# Increase recursion limit for potentially deep expression trees
sys.setrecursionlimit(20000)

class Combinator:
    """Represents a named combinator like S, K, I, or a symbolic variable."""
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return self.name

class App:
    """Represents the application of one expression to another, e.g., (f x)."""
    def __init__(self, func, arg):
        self.func = func
        self.arg = arg
    def __repr__(self):
        return f"({self.func} {self.arg})"

# Singleton combinator instances
S = Combinator("S")
K = Combinator("K")
I = Combinator("I")

def reduce_step(expr):
    """
    Performs a single, leftmost-outermost reduction step on a combinator expression.
    Returns a tuple (new_expression, was_reduced_flag).
    """
    # If not an application, it cannot be a redex (reducible expression).
    if not isinstance(expr, App):
        return expr, False

    # Try to reduce the function part of the application first (leftmost).
    new_func, reduced = reduce_step(expr.func)
    if reduced:
        return App(new_func, expr.arg), True

    # If the function part can't be reduced, try the argument part.
    new_arg, reduced = reduce_step(expr.arg)
    if reduced:
        return App(expr.func, new_arg), True

    # If no inner parts were reduced, check the current application for a redex.
    f = expr.func
    # I combinator rule: I x -> x
    if f == I:
        return expr.arg, True
    
    if isinstance(f, App):
        f_func = f.func
        # K combinator rule: K x y -> x
        if f_func == K:
            # f is (K x), expr is ((K x) y)
            x = f.arg
            return x, True
        
        if isinstance(f_func, App):
            f_func_func = f_func.func
            # S combinator rule: S x y z -> (x z) (y z)
            if f_func_func == S:
                # f_func is (S x), f is ((S x) y), expr is (((S x) y) z)
                x = f_func.arg
                y = f.arg
                z = expr.arg
                return App(App(x, z), App(y, z)), True

    # No reduction was possible at this level or below.
    return expr, False

def reduce_full(expr, max_steps=10000):
    """Fully reduces an expression by repeatedly applying reduce_step."""
    for _ in range(max_steps):
        expr, reduced = reduce_step(expr)
        if not reduced:
            break
    else:
        # This would be reached if the loop finishes without breaking.
        raise RuntimeError("Reduction exceeded maximum steps, may be a non-terminating expression.")
    return expr

def count_applications(expr, func_name):
    """Counts how many times a function `func_name` is applied in a nested expression."""
    count = 0
    current = expr
    while isinstance(current, App) and repr(current.func) == func_name:
        count += 1
        current = current.arg
    return count

def main():
    """Main function to build, evaluate the expression and find log2(n)."""
    # Build the expression: S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))
    # Standard notation `xyz` is parsed as `((x y) z)`

    # Sub-expression: P = succ = S(S(K S) K)
    P = App(S, App(App(S, App(K, S)), K))

    # Sub-expression: church_two = S(I)(I)
    church_two = App(App(S, I), I)

    # Term 2: P(church_two) = succ(2) -> should be 3
    term2 = App(P, church_two)

    # Term 1: K(P(I)) = K(succ(I)) -> K(2)
    succ_I = App(P, I)
    term1 = App(K, succ_I)

    # Let's combine them according to the full structure:
    # A(B)(C) is ((A B) C)
    # The structure is S(I) (S(I) (S(I)(Term1))(Term2))
    
    # Innermost application: S(I)(Term1)
    inner_app = App(App(S, I), term1)
    
    # Next level: (S(I)(Term1))(Term2)
    core_expr = App(inner_app, term2)
    
    # Wrapping with S(I)
    mid_expr = App(App(S, I), core_expr)
    
    # Final wrapping with S(I)
    final_expr_n = App(App(S, I), mid_expr)

    # To find the value of the Church numeral n, we evaluate n(f)(x)
    # and see how many times 'f' is applied to 'x'.
    f = Combinator("f")
    x = Combinator("x")
    
    expr_to_evaluate = App(App(final_expr_n, f), x)

    # Fully reduce the expression
    normal_form = reduce_full(expr_to_evaluate)
    
    # Count the applications of 'f'
    n = count_applications(normal_form, "f")

    # Calculate log base 2
    log2_n = math.log2(n)

    # Print the results as requested
    print(f"The Church integer n is: {n}")
    print(f"The final equation is: log_2({n}) = {log2_n}")

if __name__ == "__main__":
    main()
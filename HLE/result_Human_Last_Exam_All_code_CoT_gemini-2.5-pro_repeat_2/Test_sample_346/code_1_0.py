import math
import sys

# Set a higher recursion limit for processing deep expression trees.
sys.setrecursionlimit(2000)

# 1. Data Structures for SKI Combinators and Expressions

class Combinator:
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return self.name

S, K, I = Combinator("S"), Combinator("K"), Combinator("I")

class Var:
    def __init__(self, name):
        self.name = name
    def __repr__(self):
        return self.name

class App:
    def __init__(self, func, arg):
        self.func = func
        self.arg = arg
    def __repr__(self):
        return f"({self.func} {self.arg})"

# 2. Reduction Engine

def reduce_step(expr):
    """Performs one step of reduction using applicative order."""
    if not isinstance(expr, App):
        return expr, False

    # Try to reduce the function part of the application first
    new_func, changed_func = reduce_step(expr.func)
    if changed_func:
        return App(new_func, expr.arg), True

    # If func part didn't change, try to reduce the argument part
    new_arg, changed_arg = reduce_step(expr.arg)
    if changed_arg:
        return App(expr.func, new_arg), True
        
    # If no sub-expression changed, try to apply a rule at the top level
    f = expr.func
    x = expr.arg

    # Rule: I x -> x
    if isinstance(f, Combinator) and f.name == "I":
        return x, True

    if isinstance(f, App):
        # Rule: K x y -> x
        if isinstance(f.func, Combinator) and f.func.name == "K":
            return f.arg, True
        
        if isinstance(f.func, App):
            # Rule: S x y z -> x z (y z)
            if isinstance(f.func.func, Combinator) and f.func.func.name == "S":
                s_x = f.func.arg
                s_y = f.arg
                s_z = x # The argument from the outermost App
                return App(App(s_x, s_z), App(s_y, s_z)), True

    return expr, False # No reduction was possible

def reduce_full(expr):
    """Repeatedly applies reduction rules until in normal form."""
    while True:
        expr, changed = reduce_step(expr)
        if not changed:
            break
    return expr

# 3. Helper function to count 'f' applications in a Church numeral form
def count_applications(expr, func_var, val_var):
    """Counts n in an expression of the form f(f(...f(x)...))."""
    count = 0
    current_expr = expr
    while isinstance(current_expr, App):
        if current_expr.func == func_var:
            count += 1
            current_expr = current_expr.arg
        else:
            # Not in the expected f(...) form
            return f"Error: Unexpected function '{current_expr.func}' in result."
    
    if current_expr == val_var:
        return count
    else:
        # The chain does not end with 'x'
        return f"Error: Expression does not end with '{val_var}'."

# 4. Main script logic
if __name__ == '__main__':
    # Build the expression tree for:
    # S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))
    
    # SUCC = S(S(KS)K)
    succ_comb = App(S, App(App(S, App(K, S)), K))
    
    # A = K(SUCC(I))
    A = App(K, App(succ_comb, I))
    
    # B = SUCC(S(I)(I))
    sii = App(App(S, I), I)
    B = App(succ_comb, sii)
    
    # Expression construction from inside out
    # T = S(I)
    T = App(S, I)
    # Y = T(A) = S(I)(A)
    Y = App(T, A)
    # X = T(Y)(B) = S(I)(S(I)(A))(B)
    X = App(App(T, Y), B)
    # church_n_expr = T(X) = S(I)(S(I)(S(I)(A))(B))
    church_n_expr = App(T, X)

    # To find the value of n, we apply the numeral to symbolic 'f' and 'x'
    f = Var('f')
    x = Var('x')
    
    # This is the expression n(f)(x)
    application_expr = App(App(church_n_expr, f), x)
    
    # Reduce the full application to its normal form f(f(...x...))
    final_form = reduce_full(application_expr)
    
    # Count n from the final form
    n = count_applications(final_form, f, x)

    if isinstance(n, int):
        # Calculate log base 2 of n
        log_n = math.log2(n)
        
        print("The church integer is n.")
        print("The final equation is: log_2(n) = result")
        print("---")
        print(f"Value of n: {n}")
        print(f"Base of logarithm: 2")
        print(f"Result of log_2({n}): {log_n}")
    else:
        print(f"Could not determine n. Error: {n}")
        print(f"Final reduced form was: {final_form}")
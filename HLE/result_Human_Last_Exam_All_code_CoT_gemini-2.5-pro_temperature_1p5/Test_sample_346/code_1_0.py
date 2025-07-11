import sys
import math

# A higher recursion limit might be needed for deeply nested expressions.
sys.setrecursionlimit(20000)

# Step 1: Define classes for SKI calculus representation.
# We have a base class `Expr` and specific classes for S, K, I, and Application.

class Expr:
    """Base class for all SKI expressions."""
    def __call__(self, arg):
        # Enables a clean syntax for application, e.g., S(K)(I)
        return App(self, arg)

class S_(Expr):
    def __repr__(self): return "S"

class K_(Expr):
    def __repr__(self): return "K"

class I_(Expr):
    def __repr__(self): return "I"

class App(Expr):
    """Represents the application of one expression to another."""
    def __init__(self, func, arg):
        self.func = func
        self.arg = arg
    def __repr__(self):
        return f"({self.func} {self.arg})"

# To evaluate the final number, we need to bridge SKI logic with Python values.
class PyVal(Expr):
    """A wrapper for a Python value (e.g., an integer) inside the SKI system."""
    def __init__(self, val):
        self.val = val
    def __repr__(self):
        return f"PyVal({self.val})"

class PyFunc(Expr):
    """A wrapper for a Python function."""
    def __init__(self, func):
        self.func = func
    def __repr__(self):
        return "PyFunc"

# Create singleton instances of the basic combinators.
S, K, I = S_(), K_(), I_()

# Step 2: Build the expression tree from the problem statement.
# Expression: S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))

# For clarity, let's define the standard B (composition) and SUCC (successor) combinators.
# B = S(KS)K = λf.λg.λx. f(g(x))
B = S(K(S))(K)
# SUCC = S(B) = λn.λf.λx. f(n f x)
SUCC = S(B)

# The expression can be seen as: term_A(term_B)(term_C)
# term_A = S(I)
# term_B = S(I)(S(I)(K(SUCC(I))))
# term_C = SUCC(S(I)(I))

term_C = SUCC(S(I)(I))
term_B_inner = K(SUCC(I))
term_B_middle = S(I)(term_B_inner)
term_B = S(I)(term_B_middle)
term_A = S(I)

problem_expression = term_A(term_B)(term_C)

# Step 3: Implement the reduction logic.
def reduce_step(expr):
    """Performs a single step of normal-order (leftmost, outermost) reduction."""
    if not isinstance(expr, App):
        return expr, False # Base combinators or values cannot be reduced.

    f = expr.func
    
    # Rule: I x => x
    if isinstance(f, I_):
        return expr.arg, True
    
    # Bridge to Python: when a PyFunc is applied to a PyVal.
    if isinstance(f, PyFunc):
        if isinstance(expr.arg, PyVal):
            return PyVal(f.func(expr.arg.val)), True
        else: # Argument is not a value yet, so it must be reduced first.
            new_arg, changed = reduce_step(expr.arg)
            return (App(f, new_arg), True) if changed else (expr, False)

    if isinstance(f, App):
        ff = f.func
        # Rule: K x y => x
        if isinstance(ff, K_):
            return f.arg, True
        
        if isinstance(ff, App):
            fff = ff.func
            # Rule: S x y z => x z (y z)
            if isinstance(fff, S_):
                x, y, z = ff.arg, f.arg, expr.arg
                return x(z)(y(z)), True

    # If no top-level rule applies, try to reduce the function part of the application.
    new_func, changed = reduce_step(f)
    if changed:
        return App(new_func, expr.arg), True

    # If the function part is stable, try to reduce the argument part.
    new_arg, changed = reduce_step(expr.arg)
    if changed:
        return App(f, new_arg), True

    return expr, False # No reduction was possible.

def main():
    """Main function to evaluate the expression and find the answer."""
    # Step 4: To find the integer value of a Church numeral `n`, we compute `n(f)(x)`.
    successor_func = PyFunc(lambda i: i + 1)
    zero_val = PyVal(0)

    # The expression to evaluate is `problem_expression(f)(x)`.
    eval_expr = problem_expression(successor_func)(zero_val)
    
    # Repeatedly apply reduction until we get a final Python value.
    changed = True
    while changed and not isinstance(eval_expr, PyVal):
        eval_expr, changed = reduce_step(eval_expr)
    
    if isinstance(eval_expr, PyVal):
        n = eval_expr.val
        print(f"The SKI combinator expression evaluates to the Church integer n = {n}")
        # Step 5: Calculate and print the result.
        if n > 0:
            result = math.log2(n)
            print(f"The equation to find the answer is: log_2({n})")
            print(f"The result is: {result}")
        else:
            print(f"Cannot compute log_2({n})")
    else:
        print("The expression did not evaluate to a numerical value.")
        print(f"The final non-numeric form is: {eval_expr}")

main()
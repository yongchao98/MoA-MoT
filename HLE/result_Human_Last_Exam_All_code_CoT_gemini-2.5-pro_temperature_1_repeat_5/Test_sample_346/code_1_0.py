import sys
import math

# Set a higher recursion limit for processing deeply nested expressions.
sys.setrecursionlimit(20000)

class Combinator:
    """Base class for combinators to enable call chaining."""
    def __call__(self, arg):
        return App(self, arg)
    
    def __repr__(self):
        return self.__class__.__name__[0]

# Define the three basic combinators S, K, I.
class S_c(Combinator): pass
class K_c(Combinator): pass
class I_c(Combinator): pass

# Define symbolic placeholders for the function 'f' and value 'x'
# used to evaluate the Church numeral.
class F_c(Combinator): pass
class X_c(Combinator): pass

class App(Combinator):
    """Represents the application of one expression to another."""
    def __init__(self, p, q):
        self.p = p
        self.q = q
    
    def __repr__(self):
        return f"({self.p} {self.q})"

# Instantiate the combinators for use in parsing.
S = S_c()
K = K_c()
I = I_c()
f = F_c()
x = X_c()

def parse_expression(s: str) -> Combinator:
    """
    Parses a string representation of a combinator expression into objects.
    This uses eval() in a controlled manner, mapping names to safe objects.
    """
    # Replace combinator names with their class instantiations.
    s = s.replace('S', 'S')
    s = s.replace('K', 'K')
    s = s.replace('I', 'I')
    # Convert function call syntax `A(B)` to application syntax `A(B)`.
    s = s.replace(')(', ')(')
    
    # Evaluate the string in a controlled scope.
    return eval(s, {"S": S, "K": K, "I": I})

def reduce_step(expr: Combinator) -> (Combinator, bool):
    """
    Performs a single reduction step on an expression tree.
    Returns the new expression and a boolean indicating if a change was made.
    """
    if not isinstance(expr, App):
        return expr, False

    # Reduction Rule: I x -> x
    if isinstance(expr.p, I_c):
        return expr.q, True

    # Reduction Rule: K x y -> x
    if isinstance(expr.p, App) and isinstance(expr.p.p, K_c):
        return expr.p.q, True

    # Reduction Rule: S x y z -> x z (y z)
    if isinstance(expr.p, App) and isinstance(expr.p.p, App) and isinstance(expr.p.p.p, S_c):
        x_comb, y_comb, z_comb = expr.p.p.q, expr.p.q, expr.q
        return App(App(x_comb, z_comb), App(y_comb, z_comb)), True

    # If no top-level rule applies, reduce the left part of the application.
    new_p, changed = reduce_step(expr.p)
    if changed:
        return App(new_p, expr.q), True

    # If the left part is stable, reduce the right part.
    new_q, changed = reduce_step(expr.q)
    if changed:
        return App(expr.p, new_q), True

    # If no part can be reduced, the expression is stable for this step.
    return expr, False

def reduce_full(expr: Combinator) -> Combinator:
    """
    Repeatedly applies reduction steps until the expression is in normal form.
    """
    # A safety counter to prevent potential infinite loops in non-terminating expressions.
    max_reductions = 30000 
    for _ in range(max_reductions):
        expr, changed = reduce_step(expr)
        if not changed:
            break
    else:
        raise RuntimeError("Reduction limit exceeded, expression may not terminate.")
    return expr

def main():
    """
    Main function to execute the plan: parse, reduce, evaluate, and print the result.
    """
    expression_string = "S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))"
    
    # 1. Parse the string into a combinator object.
    church_n_combinator = parse_expression(expression_string)
    
    # 2. Apply the numeral to symbolic 'f' and 'x' to evaluate it.
    expression_to_reduce = App(App(church_n_combinator, f), x)
    
    # 3. Reduce the expression to its normal form.
    final_expression = reduce_full(expression_to_reduce)
    
    # 4. Count the applications of 'f' to determine the integer value.
    # The normal form of a Church numeral n applied to f and x is f(f(...f(x)...)).
    # We can count the occurrences of '(f ' in the string representation.
    final_expression_str = repr(final_expression)
    n_value = final_expression_str.count(f'({f} ')
    
    # 5. Calculate log_2(n) and print the results.
    log2_n = math.log2(n_value)
    
    print(f"The integer value of the Church numeral is:")
    print(f"n = {n_value}")
    print(f"The value of log_2(n) is:")
    print(f"log_2({n_value}) = {log2_n}")

if __name__ == "__main__":
    main()
import math
import sys
from functools import total_ordering

def solve():
    """
    Solves the problem by modeling SKI combinator logic in Python,
    reducing the given expression to a Church numeral, and then
    calculating the required logarithm.
    """
    # Set a higher recursion limit for deep object hashing and representation.
    sys.setrecursionlimit(20000)

    # Step 1: Define hashable data structures for combinators, variables, and applications.
    @total_ordering
    class S:
        def __repr__(self): return "S"
        def __eq__(self, other): return isinstance(other, S)
        def __lt__(self, other): return type(other).__name__ > "S"
        def __hash__(self): return hash("S")

    @total_ordering
    class K:
        def __repr__(self): return "K"
        def __eq__(self, other): return isinstance(other, K)
        def __lt__(self, other):
            if isinstance(other, S): return False
            return type(other).__name__ > "K"
        def __hash__(self): return hash("K")

    @total_ordering
    class I:
        def __repr__(self): return "I"
        def __eq__(self, other): return isinstance(other, I)
        def __lt__(self, other):
            if isinstance(other, (S, K)): return False
            return type(other).__name__ > "I"
        def __hash__(self): return hash("I")

    @total_ordering
    class Var:
        def __init__(self, name): self.name = name
        def __repr__(self): return self.name
        def __eq__(self, other): return isinstance(other, Var) and self.name == other.name
        def __lt__(self, other):
            if not isinstance(other, Var): return True
            return self.name < other.name
        def __hash__(self): return hash(self.name)

    @total_ordering
    class App:
        def __init__(self, f, x):
            self.f = f
            self.x = x
            self._hash = None
        def __repr__(self): return f"({self.f!r} {self.x!r})"
        def __eq__(self, other): return isinstance(other, App) and self.f == other.f and self.x == other.x
        def __lt__(self, other):
            if not isinstance(other, App): return True
            if self.f != other.f: return self.f < other.f
            return self.x < other.x
        def __hash__(self):
            if self._hash is None: self._hash = hash((self.f, self.x))
            return self._hash

    _reduction_cache = {}

    def reduce_expr(expr):
        """Reduces an expression to its normal form using a cache."""
        original_expr = expr
        if original_expr in _reduction_cache:
            return _reduction_cache[original_expr]

        current_expr = expr
        # Use a loop to prevent deep recursion on the top-level reductions.
        while True:
            next_expr, changed = _reduce_step(current_expr)
            if not changed:
                _reduction_cache[original_expr] = next_expr
                return next_expr
            current_expr = next_expr

    def _reduce_step(expr):
        """Performs one step of leftmost-outermost reduction."""
        if isinstance(expr, (S, K, I, Var)):
            return expr, False
        
        # Outermost reduction: Check for S, K, I redexes at the top level.
        # Rule S: S(f)(g)(x) -> f(x)(g(x))
        if isinstance(expr.f, App) and isinstance(expr.f.f, App) and isinstance(expr.f.f.f, S):
            f, g, x = expr.f.f.x, expr.f.x, expr.x
            return App(App(f, x), App(g, x)), True

        # Rule K: K(x)(y) -> x
        if isinstance(expr.f, App) and isinstance(expr.f.f, K):
            x = expr.f.x
            return x, True

        # Rule I: I(x) -> x
        if isinstance(expr.f, I):
            x = expr.x
            return x, True
        
        # If no top-level redex, reduce sub-expressions (left first for leftmost strategy).
        new_f, changed_f = _reduce_step(expr.f)
        if changed_f:
            return App(new_f, expr.x), True

        new_x, changed_x = _reduce_step(expr.x)
        if changed_x:
            return App(expr.f, new_x), True
            
        return expr, False

    def count_f_applications(expr, f_var, x_var):
        """Counts how many times `f_var` is applied in a Church numeral structure."""
        if expr == x_var:
            return 0
        if isinstance(expr, App) and expr.f == f_var:
            return 1 + count_f_applications(expr.x, f_var, x_var)
        raise ValueError(f"Expression is not a simple Church Numeral form. Final reduced form: {expr!r}")

    # Step 2: Manually and carefully construct the expression tree.
    # Define basic combinators
    s, k, i = S(), K(), I()
    
    # Build standard composite combinators from S, K, I
    # B combinator (composition): Bfgx = f(gx), B = S(KS)K
    b_comb = App(App(s, App(k, s)), k)
    # Successor function: succ nfx = f(nfx), succ = S(B)
    succ = App(s, b_comb)
    # Omega combinator (self-application): ωx = xx, ω = SII
    omega = App(App(s, i), i)
    
    # Build the main components of the given expression, assuming left-associativity
    # expr = S(I) (S(I) (S(I) (K(succ(I))) (succ(omega))))
    # The form `a b c d e` is parsed as `((((a b) c) d) e)`.
    t1 = App(s, i) # S(I)
    t2 = App(s, i) # S(I)
    t3 = App(s, i) # S(I)
    t4 = App(k, App(succ, i)) # K(succ(I))
    t5 = App(succ, omega)     # succ(omega)
    
    church_expr = App(App(App(App(t1, t2), t3), t4), t5)
    
    # Step 4: Evaluate the church numeral 'n' by applying it to symbolic f and x.
    f_var = Var('f')
    x_var = Var('x')
    
    applied_expr = App(App(church_expr, f_var), x_var)
    
    reduced_expr = reduce_expr(applied_expr)
    
    # Step 5: Count the applications to find the integer value 'n'.
    n = count_f_applications(reduced_expr, f_var, x_var)

    # Step 6: Calculate the final answer.
    result = math.log2(n)

    print(f"The integer n is {n}.")
    print(f"log_2({n}) = {result}")

solve()
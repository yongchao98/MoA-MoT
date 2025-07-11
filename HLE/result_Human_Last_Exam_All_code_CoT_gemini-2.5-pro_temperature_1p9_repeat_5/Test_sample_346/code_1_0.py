import sys
import math
import re

# Increase recursion limit for potentially deep expression trees during reduction
sys.setrecursionlimit(4000)

# --- Part 1: SKI Expression Data Structures ---

class S_class:
    """Singleton class for the S combinator."""
    def __repr__(self): return "S"

class K_class:
    """Singleton class for the K combinator."""
    def __repr__(self): return "K"

class I_class:
    """Singleton class for the I combinator."""
    def __repr__(self): return "I"

# Create singleton instances of the combinators
S, K, I = S_class(), K_class(), I_class()

class App:
    """Represents the application of a function to an argument."""
    def __init__(self, func, arg):
        self.func = func
        self.arg = arg

    def __repr__(self):
        # Use a safe representation that avoids deep recursion for unreduced trees
        return f"App({self.func!r}, {self.arg!r})"

class Val:
    """Wrapper for Python values to be used inside the SKI expression."""
    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return f"Val({self.value!r})"

# --- Part 2: Expression Parser ---

def parse_ski(s: str):
    """Parses a string like 'S(I)(K)' into a nested App object."""
    s = s.strip()
    # Base cases for atomic combinators
    if s == "S": return S
    if s == "K": return K
    if s == "I": return I
    
    # The expressions are of the form `F(A)(B)...`. We find the root function `F`
    # and then parse each parenthesized argument `(A)`, `(B)`, etc.
    try:
        first_paren_idx = s.index('(')
        func_str = s[:first_paren_idx]
        args_str = s[first_paren_idx:]
    except ValueError:
        raise ValueError(f"Invalid application format: Cannot find '(' in '{s}'")
    
    # Recursively parse the function part.
    expr = parse_ski(func_str)
    
    # Iteratively parse and apply arguments.
    # Arguments are contained in matching parentheses at the top level of args_str.
    arg_strs = []
    balance = 0
    start = 0
    if args_str:
        for i, char in enumerate(args_str):
            if char == '(':
                if balance == 0:
                    start = i + 1
                balance += 1
            elif char == ')':
                balance -= 1
                if balance == 0:
                    arg_strs.append(args_str[start:i])
    
    if balance != 0:
        raise ValueError(f"Mismatched parentheses in arguments: {args_str}")

    for arg_s in arg_strs:
        expr = App(expr, parse_ski(arg_s))
        
    return expr

# --- Part 3: SKI Expression Reducer ---

# Using a cache (memoization) to speed up reduction of common sub-expressions
reduction_cache = {}

def step(expr):
    """Performs one step of reduction using the leftmost-outermost strategy."""
    
    if not isinstance(expr, App):
        return expr # S, K, I, and Val are already in normal form
    
    expr_key = repr(expr)
    if expr_key in reduction_cache:
        return reduction_cache[expr_key]

    # Rule I: I x -> x
    if isinstance(expr.func, I_class):
        result = expr.arg
        reduction_cache[expr_key] = result
        return result

    # Rule K: (K x) y -> x
    if isinstance(expr.func, App) and isinstance(expr.func.func, K_class):
        result = expr.func.arg
        reduction_cache[expr_key] = result
        return result

    # Rule S: ((S x) y) z -> (x z) (y z)
    if isinstance(expr.func, App) and isinstance(expr.func.func, App) and isinstance(expr.func.func.func, S_class):
        x, y, z = expr.func.func.arg, expr.func.arg, expr.arg
        result = App(App(x, z), App(y, z))
        reduction_cache[expr_key] = result
        return result
        
    # Rule for applying Python callables: (Val(f)) (Val(v)) -> Val(f(v))
    if isinstance(expr.func, Val) and callable(expr.func.value) and isinstance(expr.arg, Val):
        result = Val(expr.func.value(expr.arg.value))
        reduction_cache[expr_key] = result
        return result

    # If no top-level reduction, reduce sub-expressions (leftmost-outermost)
    reduced_func = step(expr.func)
    if reduced_func is not expr.func:
        result = App(reduced_func, expr.arg)
        reduction_cache[expr_key] = result
        return result

    reduced_arg = step(expr.arg)
    if reduced_arg is not expr.arg:
        result = App(expr.func, reduced_arg)
        reduction_cache[expr_key] = result
        return result

    return expr

def reduce_full(expr):
    """Repeatedly applies the step function until a normal form is reached."""
    seen_states = set()
    while True:
        current_key = repr(expr)
        if current_key in seen_states:
             raise RuntimeError(f"Loop detected during reduction for expression: {current_key}")
        seen_states.add(current_key)

        next_expr = step(expr)
        if next_expr is expr:
            return expr # Reached normal form
        expr = next_expr

def solve():
    """
    Solves the problem by parsing and evaluating the SKI expression to find n.
    """
    expression_string = "S(I)(S(I)(S(I)(K(S(S(K(S))(K))(I)))))(S(S(K(S))(K))(S(I)(I))))"
    
    try:
        church_numeral_expr = parse_ski(expression_string)

        # To find what integer n the expression represents, we apply it to a successor
        # function (lambda x: x + 1) and a starting value (0). The result will be n.
        successor_func = Val(lambda x: x + 1)
        zero = Val(0)
        
        # Construct the test expression: ((church_numeral_expr successor_func) zero)
        test_expr = App(App(church_numeral_expr, successor_func), zero)
        
        result_val = reduce_full(test_expr)
        
        if isinstance(result_val, Val) and isinstance(result_val.value, int):
            n = result_val.value
            
            # Per instructions, outputting the numbers in the final equation.
            print(f"The complex SKI expression evaluates to the Church numeral n = {n}.")
            
            log_val = math.log2(n)

            print(f"The equation to solve is: log_2(n)")
            print(f"Substituting n = {n}, the final equation is: log_2({n}) = {log_val}")
        else:
            print(f"Error: Could not reduce the expression to an integer. Final form: {result_val!r}")

    except (ValueError, RuntimeError, RecursionError) as e:
        print(f"An error occurred: {e}")

# Execute the solution
solve()
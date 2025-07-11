def solve():
    """
    This script tests the fixed point equality fp(f*g) = fp(f) intersect fp(g)
    using a counter-example.
    The operation f*g is the pointwise meet f(x) ^ g(x).
    """

    # 1. Define the Poset (L, <=)
    # We use a 4-element Boolean lattice {0, a, b, 1}, which is distributive.
    # Elements are represented as strings for clarity.
    L = {'0', 'a', 'b', '1'}

    # Define the meet (^) and join (v) operations for L
    meet_table = {
        ('0', '0'): '0', ('0', 'a'): '0', ('0', 'b'): '0', ('0', '1'): '0',
        ('a', '0'): '0', ('a', 'a'): 'a', ('a', 'b'): '0', ('a', '1'): 'a',
        ('b', '0'): '0', ('b', 'a'): '0', ('b', 'b'): 'b', ('b', '1'): 'b',
        ('1', '0'): '0', ('1', 'a'): 'a', ('1', 'b'): 'b', ('1', '1'): '1',
    }
    join_table = {
        ('0', '0'): '0', ('0', 'a'): 'a', ('0', 'b'): 'b', ('0', '1'): '1',
        ('a', '0'): 'a', ('a', 'a'): 'a', ('a', 'b'): '1', ('a', '1'): '1',
        ('b', '0'): 'b', ('b', 'a'): '1', ('b', 'b'): 'b', ('b', '1'): '1',
        ('1', '0'): '1', ('1', 'a'): '1', ('1', 'b'): '1', ('1', '1'): '1',
    }
    
    def meet(x, y):
        return meet_table[(x, y)]

    def join(x, y):
        return join_table[(x, y)]

    # 2. Define functions f and g
    # f(x) = x v a
    # g(x) = x v b
    # These functions are extensive (x <= f(x), x <= g(x)), monotone, and continuous.
    def f(x):
        return join(x, 'a')

    def g(x):
        return join(x, 'b')

    # 3. Compute the set of fixed points for f and g
    fp_f = {x for x in L if f(x) == x}
    fp_g = {x for x in L if g(x) == x}

    # 4. Compute the right-hand side of the equation: fp(f) intersect fp(g)
    RHS = fp_f.intersection(fp_g)

    # 5. Define the function f*g (pointwise meet)
    def f_meet_g(x):
        return meet(f(x), g(x))

    # 6. Compute the left-hand side of the equation: fp(f*g)
    LHS = {x for x in L if f_meet_g(x) == x}
    
    # 7. Print the results
    print("--- Analysis of fp(f*g) = fp(f) intersect fp(g) ---")
    print(f"Poset L = {sorted(list(L))}")
    print("f(x) = x v 'a', g(x) = x v 'b'")
    print("f and g are extensive, monotone, and continuous.")
    print("")

    # Print the equation components
    print(f"fp(f) = {sorted(list(fp_f))}")
    print(f"fp(g) = {sorted(list(fp_g))}")
    print(f"The set of fixed points of f*g is: fp(f*g) = {sorted(list(LHS))}")
    print(f"The intersection of fixed points is: fp(f) ∩ fp(g) = {sorted(list(RHS))}")

    # Final conclusion from the test
    print("\nFinal Equation Check:")
    print(f"{sorted(list(LHS))} = {sorted(list(RHS))}")
    if LHS == RHS:
        print("The equality holds for this example.")
    else:
        print("The equality does NOT hold for this example.")
        print("This disproves that 'f and g extensive/monotone/continuous' is a sufficient condition.")

solve()

# Note on a possible interpretation: If f*g were meant to be the pointwise JOIN (f V g),
# the minimal requirement for the equality fp(f V g) = fp(f) intersect fp(g) to hold
# is that both f and g are extensive. However, the dot notation strongly suggests meet.

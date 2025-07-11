import itertools

def solve():
    """
    Calculates the number of unique light ray paths by determining the number of
    unique final virtual images formed by the sequence of reflections.
    """

    # A function to compose two linear transformations of the form f(z) = az + b.
    # A transformation is represented by a tuple (a, b).
    # If we have g(f(z)), where f(z) = a*z + b and g(z') = c*z' + d,
    # then g(f(z)) = c*(a*z + b) + d = (c*a)*z + (c*b + d).
    # The composition of g after f is represented by the tuple (c*a, c*b + d).
    def compose(g_op, f_op):
        c, d = g_op
        a, b = f_op
        return (c * a, c * b + d)

    print("Step 1: Analyze horizontal reflections (2 on G1, 1 on G3).")
    # T1 (G1): x -> -x. Represented as (-1, 0).
    # T3 (G3): x -> 2w - x. Represented as (-1, 2) (using a normalized width w=1).
    T1 = (-1, 0)
    T3 = (-1, 2)
    horizontal_ops = [T1, T1, T3]
    
    # Get all unique sequences (permutations) of the horizontal operators.
    unique_h_permutations = sorted(list(set(itertools.permutations(horizontal_ops))))
    
    horizontal_outcomes = set()
    for p in unique_h_permutations:
        # The path involves reflections in sequence, e.g., op1 then op2 then op3.
        # The final virtual image is the result of op3(op2(op1(M))).
        # We start with the identity transformation, z -> 1*z + 0, represented by (1, 0).
        final_op = (1, 0)
        for op in p:
            final_op = compose(op, final_op)
        horizontal_outcomes.add(final_op)

    num_h_outcomes = len(horizontal_outcomes)
    print(f"The possible sequences of horizontal reflections produce {num_h_outcomes} unique outcomes.")

    print("\nStep 2: Analyze vertical reflections (1 on G2, 1 on G4).")
    # T2 (G2): y -> -y. Represented as (-1, 0).
    # T4 (G4): y -> 2h - y. Represented as (-1, 2) (using a normalized height h=1).
    T2 = (-1, 0)
    T4 = (-1, 2)
    vertical_ops = [T2, T4]

    # Get all unique sequences (permutations) of the vertical operators.
    unique_v_permutations = sorted(list(set(itertools.permutations(vertical_ops))))
    
    vertical_outcomes = set()
    for p in unique_v_permutations:
        final_op = (1, 0)
        for op in p:
            final_op = compose(op, final_op)
        vertical_outcomes.add(final_op)

    num_v_outcomes = len(vertical_outcomes)
    print(f"The possible sequences of vertical reflections produce {num_v_outcomes} unique outcomes.")
    
    print("\nStep 3: Calculate the total number of unique paths.")
    print("The total number of paths is the product of the number of unique horizontal and vertical outcomes.")
    
    total_ways = num_h_outcomes * num_v_outcomes
    
    print(f"Number of ways = (Unique horizontal outcomes) * (Unique vertical outcomes)")
    print(f"Number of ways = {num_h_outcomes} * {num_v_outcomes} = {total_ways}")

# Run the solver
solve()
print("<<<4>>>")
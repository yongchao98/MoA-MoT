import sys

def solve_torus_action():
    """
    This script demonstrates the two actions of the fundamental group of a torus
    on the fiber of its universal cover and shows they are equivalent.
    """
    
    # 1. Model the mathematical objects for the torus T^2.
    # The fundamental group pi_1(T^2) is Z^2.
    # The fiber p^-1(x_0) is also Z^2.
    
    print("Let X be the torus T^2.")
    print("The fundamental group pi_1(T^2) is isomorphic to Z^2.")
    print("The fiber of the universal cover p^-1(x_0) is also isomorphic to Z^2.")
    print("-" * 40)
    
    # 2. Choose specific examples for a loop class g and a fiber point f.
    g_loop_class = (2, 3)  # Represents an element [gamma] in pi_1(T^2)
    f_fiber_point = (4, 5) # Represents a point ~x in the fiber p^-1(x_0)
    
    m, n = g_loop_class
    k, l = f_fiber_point
    
    print(f"We choose a loop class g = ({m}, {n}) from the fundamental group.")
    print(f"We choose a point f = ({k}, {l}) in the fiber.")
    print("-" * 40)
    
    # 3. Simulate Action 1 (Holonomy / Path Lifting)
    print("--- Action 1: Holonomy around loops ---")
    print("This action is defined by lifting the loop g starting from the fiber point f.")
    print("The result is the endpoint of this lifted path.")
    print("For g=(m,n) and f=(k,l), the lifted path starts at (k,l) and ends at (k+m, l+n).")
    
    res1_k = k + m
    res1_l = l + n
    result1 = (res1_k, res1_l)
    
    print("\nCalculation:")
    print(f"The first component is {k} + {m} = {res1_k}")
    print(f"The second component is {l} + {n} = {res1_l}")
    print(f"Result of Action 1: ({res1_k}, {res1_l})")
    print("-" * 40)

    # 4. Simulate Action 2 (Deck Transformations)
    print("--- Action 2: Restricting deck transformations ---")
    print("This action is defined by applying the deck transformation phi_g, which corresponds to g, to the point f.")
    # The deck transformation phi_g is the one that maps the base fiber point (0,0) to g.
    # For the torus, this is simply the translation by the vector g.
    print("The deck transformation phi_g for g=(m,n) is the translation by the vector (m,n).")
    print("We apply this translation to the point f=(k,l).")
    
    res2_k = k + m
    res2_l = l + n
    result2 = (res2_k, res2_l)
    
    print("\nCalculation:")
    print(f"The first component is {k} + {m} = {res2_k}")
    print(f"The second component is {l} + {n} = {res2_l}")
    print(f"Result of Action 2: ({res2_k}, {res2_l})")
    print("-" * 40)

    # 5. Compare and Conclude
    print("--- Conclusion ---")
    if result1 == result2:
        print(f"Both actions yield the same result: {result1}.")
        print("This is because the fundamental group of the torus, Z^2, is abelian.")
        print("For an abelian group, the standard left and right actions are equivalent.")
        print("\nTherefore, the two actions are the same when X = T^2.")
    else:
        # This case won't be reached for the torus.
        print("The actions yield different results.")
        print(f"Action 1 Result: {result1}")
        print(f"Action 2 Result: {result2}")
        
if __name__ == "__main__":
    solve_torus_action()
    # Finally, printing the answer in the required format
    # The original question is a Yes/No question.
    sys.stdout.write("\n<<<Yes>>>\n")

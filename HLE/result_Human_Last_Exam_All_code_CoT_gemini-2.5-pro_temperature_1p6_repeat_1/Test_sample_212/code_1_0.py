import numpy as np

def holonomy_action(loop_element, fiber_point):
    """
    Calculates the action of a loop on a fiber point via holonomy (path lifting).
    For the torus, this corresponds to adding the loop's integer vector representation.
    """
    m, n = loop_element
    k, l = fiber_point
    
    # The lift of the loop (m,n) starting at (k,l) is the path t -> (mt+k, nt+l).
    # The endpoint at t=1 is (m+k, n+l).
    result = (k + m, l + n)
    
    print("1. Holonomy Action:")
    print(f"   Loop (m,n) = {loop_element}")
    print(f"   Fiber point (k,l) = {fiber_point}")
    print(f"   The action on the fiber point is vector addition.")
    print(f"   Result: {fiber_point} + {loop_element} = ({k}+{m}, {l}+{n}) = {result}")
    return result

def deck_transformation_action(loop_element, fiber_point):
    """
    Calculates the action of a loop on a fiber point via deck transformations.
    """
    m, n = loop_element
    k, l = fiber_point
    
    # Step 1: Find the deck transformation corresponding to the loop (m,n).
    # This is the deck transformation that maps the basepoint of the fiber (0,0)
    # to the endpoint of the lift of the loop from (0,0).
    # The endpoint of the lift is (m,n).
    # Deck transformations are g_ab(u,v) = (u+a, v+b).
    # We need g_ab(0,0) = (m,n), which implies (a,b) = (m,n).
    # So the corresponding deck transformation is translation by (m,n).
    deck_translation_vector = loop_element
    
    # Step 2: Apply this deck transformation to the fiber point (k,l).
    result = (k + deck_translation_vector[0], l + deck_translation_vector[1])
    
    print("\n2. Deck Transformation Action:")
    print(f"   Loop (m,n) = {loop_element}")
    print(f"   Fiber point (k,l) = {fiber_point}")
    print(f"   The corresponding deck transformation is translation by {deck_translation_vector}.")
    print(f"   Applying this transformation: {fiber_point} + {deck_translation_vector} = ({k}+{m}, {l}+{n}) = {result}")
    return result

if __name__ == '__main__':
    # Define an example loop from the fundamental group Z^2
    # This loop goes 3 times longitudinally and 5 times latitudinally.
    example_loop = (3, 5)
    
    # Define an example point in the fiber Z^2
    example_fiber_point = (1, 2)
    
    # Calculate the result of each action
    holonomy_result = holonomy_action(example_loop, example_fiber_point)
    deck_result = deck_transformation_action(example_loop, example_fiber_point)
    
    # Compare the results
    print("\n--- Comparison ---")
    if holonomy_result == deck_result:
        print(f"The results are identical: {holonomy_result}")
        print("Conclusion: For the 2-torus, these two actions are the same.")
    else:
        print(f"The results are different: Holonomy={holonomy_result}, Deck={deck_result}")
        print("Conclusion: For the 2-torus, these two actions are NOT the same.")

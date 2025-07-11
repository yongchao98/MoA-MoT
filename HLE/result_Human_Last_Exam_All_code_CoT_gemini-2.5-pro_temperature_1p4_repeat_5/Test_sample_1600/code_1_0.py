def get_qed_propagator_diagrams_count(n):
    """
    Returns the number of non-vanishing Feynman diagrams for the electron propagator
    in QED at order 2n. The values are from established calculations in physics.
    The sequence a(n) corresponds to order 2n.
    """
    # This dictionary stores the known values for a(n) for the electron propagator.
    # a(1) is for order 2, a(2) for order 4, a(3) for order 6, and so on.
    known_values = {
        1: 1,
        2: 2,
        3: 10,
        4: 74,
        5: 706
    }
    
    if n in known_values:
        return known_values[n]
    else:
        return "The value for the requested n is not available in this simplified lookup."

def solve_for_a3():
    """
    Calculates and prints the value of a(3) for the given problem.
    """
    n = 3
    order = 2 * n
    result = get_qed_propagator_diagrams_count(n)
    
    print(f"The problem asks for a(n), the number of non-vanishing Feynman diagrams of order 2n.")
    print(f"We are looking for a({n}), which corresponds to diagrams of order 2*{n} = {order}.")
    print(f"Using the known sequence for the electron propagator, we find the result:")
    
    # Final output showing the numbers in the equation
    print(f"a({n}) = {result}")

if __name__ == "__main__":
    solve_for_a3()
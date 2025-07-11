def solve_problem():
    """
    This function calculates the value of the expression sum(k * n_k for k in 1..6)
    based on the derived exponents for the magnetic field noise problem.
    """
    
    # Exponents derived from the physics of the problem
    n1 = 1
    n2 = 1
    n3 = -4
    n4 = 0
    n5 = -2
    n6 = -1
    
    exponents = [n1, n2, n3, n4, n5, n6]
    
    # Calculate the sum
    total = 0
    calculation_steps = []
    for k, n in enumerate(exponents, 1):
        term = k * n
        total += term
        
        # Format the term for the equation string
        if n < 0:
            calculation_steps.append(f"{k}*({n})")
        else:
            calculation_steps.append(f"{k}*{n}")
            
    # Build the full equation string
    equation_str = " + ".join(calculation_steps)
    
    # Print the equation and the final result
    print("The final calculation is:")
    print(f"{equation_str} = {total}")

solve_problem()
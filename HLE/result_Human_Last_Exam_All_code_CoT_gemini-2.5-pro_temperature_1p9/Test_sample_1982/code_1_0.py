def solve_control_problem():
    """
    Solves for the control variable u1 based on the given matrix equation and parameters.
    
    The plan is as follows:
    1.  The matrix equation is:
        [[0, 1], [0, 0]] * [[x11, x12], [x21, x22]] + [[x11, x12], [x21, x22]] * [[1, 0], [0, 1]] 
        = [[1, 0], [0, 1]] + [[c1, 0], [0, c2]] * [[u1, 0], [0, u2]]

    2.  Performing the matrix multiplication and addition on both sides gives:
        LHS = [[x21 + x11, x22 + x12], [x21, x22]]
        RHS = [[1 + c1*u1, 0], [0, 1 + c2*u2]]

    3.  Equating the elements of the matrices gives a system of four equations.
        From the top-left element, we get: x11 + x21 = 1 + c1*u1
        From the bottom-left element, we get: x21 = 0
        Substituting x21 = 0 into the first equation yields: x11 = 1 + c1*u1

    4.  The problem provides values for l1 and alpha1. We assume l1 corresponds to x11.
        So, we have the equation: l1 = 1 + c1 * u1.

    5.  Solving for u1, we get: u1 = (l1 - 1) / c1.
        Now we can substitute the given values to find u1.
    """
    
    # Given values
    # c1 = 10^4
    c1 = 10**4
    
    # l1 = (1 + 10^5)^5
    l1 = (1 + 10**5)**5
    
    # Calculate u1 using the formula u1 = (l1 - 1) / c1
    # Python's integers can handle the large numbers involved.
    # We use integer division // as the result is expected to be an integer.
    u1 = (l1 - 1) // c1
    
    # Output the final equation with all numbers filled in
    print("The final equation is l1 = 1 + c1 * u1")
    print("Substituting the given and calculated values:")
    # We use f-strings for clear and easy formatting
    print(f"{l1} = 1 + {c1} * {u1}")
    
    # Print the final answer for u1 as requested.
    print("\nThe calculated value for the control u1 is:")
    print(u1)
    
    # Store the final answer in the specified format for extraction
    final_answer = u1
    return final_answer

# Run the function to get the solution
final_answer = solve_control_problem()
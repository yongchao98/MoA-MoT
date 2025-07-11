import math

def get_blowup_condition(x0):
    """
    For the given system of ODEs and an initial value x(0)=x0,
    this function calculates and prints the condition on y(0)
    for the solution to blow up.
    
    The final equation is the inequality y(0) < threshold. This function
    prints the derivation of the threshold with all numerical values.
    """
    print(f"Analyzing blow-up condition for the initial value x(0) = {x0}\n")
    
    if x0 <= 1:
        print("This analysis is only valid for x(0) > 1.")
        return

    # The threshold for blow-up is determined by the separatrix connected to the saddle point (1,0).
    # The equation for this separatrix is y^2 = 2*x + 1 - 3*x^(2/3).
    # Blow-up occurs if y(0) is below the threshold y_threshold = sqrt(2*x(0) + 1 - 3*x(0)^(2/3)).

    # Calculate the components for the threshold formula.
    
    # Step 1: Calculate x0^(2/3)
    x0_pow_2_3 = x0**(2.0/3.0)
    
    # Step 2: Calculate 2*x0
    term_2x0 = 2 * x0
    
    # Step 3: Calculate 3 * x0^(2/3)
    term_3_pow = 3 * x0_pow_2_3
    
    # Step 4: Calculate the expression inside the square root
    val_in_sqrt = term_2x0 + 1 - term_3_pow
    
    # Final threshold value
    y_threshold = math.sqrt(val_in_sqrt)

    print("The condition for the solution to blow up is y(0) < Threshold.")
    print("The threshold is calculated from the separatrix equation: y^2 = 2*x + 1 - 3*x^(2/3).\n")
    
    print("Calculation of the threshold value:")
    # Printing the equation step-by-step with each number.
    print(f"Threshold = sqrt( 2 * {x0} + 1 - 3 * ({x0})^(2/3) )")
    print(f"            = sqrt( {term_2x0:.4f} + 1 - 3 * {x0_pow_2_3:.4f} )")
    print(f"            = sqrt( {term_2x0:.4f} + 1 - {term_3_pow:.4f} )")
    print(f"            = sqrt( {term_2x0 + 1:.4f} - {term_3_pow:.4f} )")
    print(f"            = sqrt( {val_in_sqrt:.4f} )")
    print(f"            = {y_threshold:.4f}\n")
    
    print("Final Answer: The solution of the system blows up if")
    print(f"y(0) < {y_threshold:.4f}")


# As the problem states x(0) > 1 but does not provide a specific value,
# we will use an example value of x(0) = 8 to demonstrate the code.
# You can change this value to any number greater than 1.
initial_x = 8.0
get_blowup_condition(initial_x)
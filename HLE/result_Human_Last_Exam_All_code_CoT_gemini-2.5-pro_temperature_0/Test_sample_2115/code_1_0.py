import numpy as np

def solve_integral():
    """
    Calculates the value of the integral based on the analytical solution.
    The analytical solution is 3 * ln(3 / (e^2 + e + 1)).
    This function breaks down the calculation and prints each step.
    """
    
    # The final analytical equation is Result = C1 * ln(C2 / (e^C3 + e + C4))
    # Here, C1=3, C2=3, C3=2, C4=1.
    c1 = 3.0
    c2 = 3.0
    c3 = 2.0
    c4 = 1.0
    
    # Value of the mathematical constant e
    e = np.e
    
    # Step-by-step calculation
    e_squared = np.power(e, c3)
    denominator = e_squared + e + c4
    fraction = c2 / denominator
    log_of_fraction = np.log(fraction)
    final_result = c1 * log_of_fraction
    
    # Output the breakdown of the final equation as requested
    print("The final analytical solution is of the form: c1 * ln(c2 / (e^c3 + e + c4))")
    print(f"The numbers in the equation are:")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    print("\n--- Calculation Steps ---")
    print(f"e = {e}")
    print(f"e^{int(c3)} = {e_squared}")
    print(f"Denominator (e^{int(c3)} + e + {int(c4)}) = {denominator}")
    print(f"Fraction ({int(c2)} / Denominator) = {fraction}")
    print(f"Natural Log of Fraction (ln({fraction:.4f})) = {log_of_fraction}")
    print(f"Final Result ({int(c1)} * {log_of_fraction:.4f}) = {final_result}")
    
    print("\n--- Final Answer ---")
    print(final_result)

solve_integral()
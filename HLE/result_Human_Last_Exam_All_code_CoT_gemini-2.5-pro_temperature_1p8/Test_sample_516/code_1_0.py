import sympy

def calculate_expected_t():
    """
    Calculates the expected value E[T] by defining and evaluating
    the symbolic infinite sum derived from the problem description.
    """
    # Define the symbolic variable for summation
    i = sympy.Symbol('i', integer=True, positive=True)

    # The general term for the sum is (1/3^i) * E[Y_i * I(m>=i)].
    # We derived that E[Y_i * I(m>=i)] = i / (i! * (i+1)).
    # We can rewrite this term as 1/((i-1)! * (i+1)) for the sum starting at i=1.
    term = (sympy.S(1)/3)**i / (sympy.factorial(i-1) * (i+1))

    # Define the infinite sum for E[T]
    sum_expr = sympy.Sum(term, (i, 1, sympy.oo))

    # Calculate the sum symbolically
    result = sum_expr.doit()

    # The final equation is E[T] = 3 - 2*exp(1/3) which is in the form A + B*exp(C/D)
    A = 3
    B = -2
    C = 1
    D = 3

    # Print the final result and the numbers in the equation
    print(f"The symbolic expression for the expected value E[T] is: {result}")
    print("\nThe final equation is of the form: A + B * exp(C/D)")
    print("The numbers in this final equation are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    
    # Print the final equation with the numbers plugged in
    print(f"\nThus, the final equation is: E[T] = {A} + ({B}) * exp({C}/{D})")

if __name__ == '__main__':
    calculate_expected_t()

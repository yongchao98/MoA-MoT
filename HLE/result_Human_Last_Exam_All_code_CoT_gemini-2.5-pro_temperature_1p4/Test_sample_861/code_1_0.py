import numpy as np

def solve_difference_equation():
    """
    Solves the given difference equation to find the constants and calculates the final expression.
    """
    # Step 1: Find the particular solution constant E
    # The equation is 8y[n] - 6y[n-1] + y[n-2] = 1.
    # For a constant particular solution y_p[n] = E, we have:
    # 8*E - 6*E + E = 1
    # 3*E = 1
    E = 1/3

    # Step 2: Find the roots of the characteristic equation for the homogeneous part
    # The characteristic equation is 8*r^2 - 6*r + 1 = 0.
    # We use the quadratic formula r = (-b +/- sqrt(b^2 - 4ac)) / (2a)
    a, b, c = 8, -6, 1
    delta = np.sqrt(b**2 - 4*a*c)
    r1 = (-b + delta) / (2*a)
    r2 = (-b - delta) / (2*a)

    # Assign roots to B and D as per the problem's format y[n] = A*(B)^n + C*(D)^n + E
    # We can assign the larger root to B.
    B = max(r1, r2)
    D = min(r1, r2)

    # Step 3: Use initial conditions to find A and C
    # The general solution is y[n] = A*B^n + C*D^n + E
    # Initial conditions: y[0] = 1, y[-1] = 2

    # From y[0] = 1:
    # A*B^0 + C*D^0 + E = 1  => A + C = 1 - E
    
    # From y[-1] = 2:
    # A*B^-1 + C*D^-1 + E = 2 => A/B + C/D = 2 - E

    # We have a system of two linear equations:
    # 1*A + 1*C = 1 - E
    # (1/B)*A + (1/D)*C = 2 - E
    
    # We can solve this using numpy.linalg.solve(M, v) for M*x = v
    # where x = [A, C]
    M = np.array([[1, 1], 
                  [1/B, 1/D]])
    v = np.array([1 - E, 
                  2 - E])
    
    solution = np.linalg.solve(M, v)
    A = solution[0]
    C = solution[1]

    print("The final closed-form solution is y[n] = A*(B)^n + C*(D)^n + E")
    print(f"The constants are:")
    print(f"A = {A}")
    print(f"B = {B}")
    print(f"C = {C}")
    print(f"D = {D}")
    print(f"E = {E}")
    print("-" * 20)

    # Step 4: Calculate the final expression E/A + (D*C)/B
    term1 = E / A
    term2 = (D * C) / B
    result = term1 + term2
    
    print("Calculating the final expression: E/A + (D*C)/B")
    print(f"E/A = ({E}) / ({A}) = {term1}")
    print(f"(D*C)/B = (({D}) * ({C})) / ({B}) = {term2}")
    print(f"Final Result = {term1} + {term2} = {result}")

solve_difference_equation()
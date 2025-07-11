import numpy as np

def solve():
    """
    Calculates the minimal polynomial of the connective constant for the given graph.

    The method involves using a transfer matrix representation of the graph. The
    characteristic polynomial for the eigenvalues (mu) of the transfer matrix is
    dependent on a wave number k:
      mu^2 - (4*cos(k))*mu + (4*cos(k)^2 - 2 - 2*cos(k)) = 0

    For the specific geometry of this graph (a strip of the triangular lattice),
    the relevant value of k is pi/3. Substituting this value gives the
    minimal polynomial for the connective constant.
    """
    # The special value of the wave number for this lattice
    k = np.pi / 3
    
    cos_k = np.cos(k)  # This is exactly 0.5
    
    # The polynomial is of the form a*x^2 + b*x + c = 0
    # The coefficient of x^2 is 1.
    a = 1.0
    
    # The coefficient of x is -4*cos(k)
    b = -4 * cos_k
    
    # The constant term is 4*cos(k)^2 - 2 - 2*cos(k)
    c = 4 * cos_k**2 - 2 - 2 * cos_k
    
    # The coefficients are integers
    a_int = int(round(a))
    b_int = int(round(b))
    c_int = int(round(c))

    # Format the output string for the polynomial equation
    # Handle signs and coefficients of 1 for nice printing
    poly_string = ""
    if a_int == 1:
        poly_string += "x^2 "
    elif a_int == -1:
        poly_string += "-x^2 "
    else:
        poly_string += f"{a_int}x^2 "

    if b_int > 0:
        poly_string += f"+ {b_int}x "
    elif b_int < 0:
        poly_string += f"- {abs(b_int)}x "

    if c_int > 0:
        poly_string += f"+ {c_int} "
    elif c_int < 0:
        poly_string += f"- {abs(c_int)} "
        
    poly_string += "= 0"

    print("The minimal polynomial for the connective constant is:")
    print(poly_string)
    
    print("\nThe coefficients of the polynomial are:")
    print(f"Coefficient of x^2: {a_int}")
    print(f"Coefficient of x: {b_int}")
    print(f"Constant term: {c_int}")

solve()
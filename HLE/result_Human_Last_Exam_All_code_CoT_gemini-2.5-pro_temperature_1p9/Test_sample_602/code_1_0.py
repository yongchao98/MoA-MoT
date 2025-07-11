import sys
import math

def calculate_l(n: int):
    """
    This function calculates the value of l(n) based on the derived analytical formula.
    
    The derivation proceeds as follows:
    1. Define the matrices M(n) and P(n).
       - M(n) is a matrix with a on the diagonal and b on the off-diagonal, where 
         a = sqrt(1 - (n-1)/n^2) and b = 1/n.
       - P(n) is simplified to be D * T_inv * D where D is a diagonal matrix of (-1)^i 
         and T is a tridiagonal matrix with 2s on diagonal and -1s on off-diagonals.
    2.  Calculate f_3(P(n)): This is the matrix inverse of P(n), which results in a matrix Y
        with 2s on the diagonal and 1s on the super- and sub-diagonals.
    3.  Calculate f_2(Y): This is the projection of Y onto the tangent space at M(n). The formula for the 
        projected matrix Z is Y - M * diag(M^T * Y).
    4.  Calculate f_1(Z): This is the sum of the elements in the first and last rows of Z.
    5.  The final exact expression for l(n) after algebraic simplification is:
        l(n) = (2*n^2 + 2 - (4*n - 2) * sqrt(n^2 - n + 1)) / n^2

    The code implements this final formula.
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    # Components of the derived formula l(n) = (A - B * sqrt(C)) / D
    A = 2 * n**2 + 2
    B = 4 * n - 2
    C = n**2 - n + 1
    D = n**2
    
    print(f"For n = {n}, we calculate l(n) using the formula: (A - B * sqrt(C)) / D")
    print("Where the components are:")
    print(f"  A = 2*n^2 + 2 = {A}")
    print(f"  B = 4*n - 2 = {B}")
    print(f"  C = n^2 - n + 1 = {C}")
    print(f"  D = n^2 = {D}")
    
    sqrt_C = math.sqrt(C)
    
    result = (A - B * sqrt_C) / D
    
    print(f"\nThe equation is: ({A} - {B} * sqrt({C})) / {D}")
    print(f"Which evaluates to: ({A} - {B} * {sqrt_C}) / {D}")
    print("\nFinal Result:")
    print(f"l({n}) = {result}")
    return result

if __name__ == '__main__':
    if len(sys.argv) > 1:
        try:
            n_val = int(sys.argv[1])
        except ValueError:
            print("Error: Please provide a valid integer for n.")
            sys.exit(1)
    else:
        # Default to n=5 if no argument is given, as n must be >= 5.
        n_val = 5
        print(f"No value for n provided. Using default value n = 5.")
        
    calculate_l(n_val)

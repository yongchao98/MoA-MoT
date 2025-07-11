import math
import numpy as np

def get_continued_fraction_coeffs(p, q):
    """
    Computes the partial quotients of the continued fraction of p/q.
    p/q is a positive rational number less than 1.
    """
    coeffs = []
    # Using the Euclidean algorithm to find the quotients
    # for q/p, which will be [a1, a2, ..., an]
    num, den = q, p
    while den != 0:
        quotient = num // den
        remainder = num % den
        coeffs.append(quotient)
        num = den
        den = remainder
    return coeffs

def format_equation(a, b, c):
    """
    Formats the coefficients a, b, c into a standard quadratic equation string.
    e.g., Ax^2 + Bx + C = 0
    """
    parts = []
    # Term Ax^2
    if a != 0:
        if a == 1:
            parts.append("x^2")
        elif a == -1:
            parts.append("-x^2")
        else:
            parts.append(f"{a}x^2")
    
    # Term Bx
    if b != 0:
        # Add a sign if it's not the first term
        if parts:
            if b > 0:
                parts.append(" + ")
            else:
                parts.append(" - ")
        # If it is the first term and negative
        elif b < 0:
            parts.append("-")
            
        b_abs = abs(b)
        if b_abs == 1:
            parts.append("x")
        else:
            parts.append(f"{b_abs}x")

    # Term C
    if c != 0:
        if parts:
            if c > 0:
                parts.append(" + ")
            else:
                parts.append(" - ")
        elif c < 0:
             parts.append("-")

        c_abs = abs(c)
        parts.append(f"{c_abs}")

    parts.append(" = 0")
    return "".join(parts)
    
def main():
    """
    Main function to solve the problem.
    """
    # Step 1: Define the rational number
    p, q = 4, 7

    # Step 2: Compute the continued fraction of p/q
    # 4/7 = [0; 1, 1, 3]
    # The partial quotients are [1, 1, 3]
    coeffs = get_continued_fraction_coeffs(p, q)
    
    # Step 3: Reverse the sequence of quotients
    # The period of the associated continued fraction is the reverse
    reversed_coeffs = list(reversed(coeffs))
    
    # Step 4: Compute the matrix K for the reversed sequence
    # K = M_{a_n} * ... * M_{a_1}
    # For period [3, 1, 1]
    K = np.identity(2, dtype=object) # Use object to handle large integers if necessary
    
    # Initialize with the first matrix in the sequence
    if reversed_coeffs:
        a = reversed_coeffs[0]
        K = np.array([[a, 1], [1, 0]], dtype=object)
        
        # Multiply by the rest
        for i in range(1, len(reversed_coeffs)):
            a = reversed_coeffs[i]
            M_a = np.array([[a, 1], [1, 0]], dtype=object)
            K = np.dot(K, M_a)

    # Step 5: Derive the quadratic equation coefficients
    # K = [[P, S], [Q, R]]
    # The equation is Qx^2 + (R-P)x - S = 0
    P, S = K[0, 0], K[0, 1]
    Q, R = K[1, 0], K[1, 1]
    
    A = int(Q)
    B = int(R - P)
    C = int(-S)
    
    # Step 6: Simplify the coefficients
    common_divisor = math.gcd(math.gcd(A, B), C)
    A //= common_divisor
    B //= common_divisor
    C //= common_divisor
    
    # Step 7: Print the final equation, which represents the continued fraction
    equation_str = format_equation(A, B, C)
    
    print(f"The continued fraction of {p}/{q} gives the quotients: {coeffs}")
    print(f"The associated periodic continued fraction has the reversed period: {reversed_coeffs}")
    print("The quadratic irrational number x with this continued fraction satisfies the equation:")
    print(equation_str)
    print(f"The numbers in the final equation are: A={A}, B={B}, C={C}")

if __name__ == "__main__":
    main()
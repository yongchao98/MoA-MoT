import sys

# It is recommended to use python3 for this script.
if sys.version_info[0] < 3:
    sys.exit("This script requires Python 3.")

def solve():
    """
    This function calculates the value of a_{n,k,l} mod p based on the given recurrence and parameters.
    """
    p = 21023

    # Coefficients from P(x,y) = 12 + 3x + 75y + 27x^2y^2
    A, B, C, D = 12, 3, 75, 27

    print(f"The prime modulus is p = {p}.")
    print("The problem reduces to calculating C^E mod p, where C is a product of terms and E is an exponent.")
    print("C = a_{5,2,2} * a_{3,1,2} * a_{2,1,1} mod p")
    print(f"E = (3*p + 1) / 2 = (3*{p} + 1) / 2")
    print("-" * 30)
    print("Calculating the components of C:")
    
    # Calculate a_{2,1,1} = [x^1 y^1] P(x,y)^2 = 2*B*C mod p
    a_2_1_1 = (2 * B * C) % p
    print(f"a_{{2,1,1}} = (2 * {B} * {C}) mod {p} = {a_2_1_1}")

    # Calculate a_{3,1,2} = [x^1 y^2] P(x,y)^3 = 3*B*C^2 mod p
    a_3_1_2 = (3 * B * pow(C, 2, p)) % p
    print(f"a_{{3,1,2}} = (3 * {B} * {C}^2) mod {p} = {a_3_1_2}")

    # Calculate a_{5,2,2} = [x^2 y^2] P(x,y)^5 = (5*A^4*D + 30*A*B^2*C^2) mod p
    term1 = (5 * pow(A, 4, p) * D) % p
    term2 = (30 * A * pow(B, 2, p) * pow(C, 2, p)) % p
    a_5_2_2 = (term1 + term2) % p
    print(f"a_{{5,2,2}} = (5*{A}^4*{D} + 30*{A}*{B}^2*{C}^2) mod {p} = {a_5_2_2}")
    
    print("-" * 30)

    # Calculate C = a_{2,1,1} * a_{3,1,2} * a_{5,2,2} mod p
    C_prod = (a_2_1_1 * a_3_1_2 * a_5_2_2) % p
    print(f"C = ({a_2_1_1} * {a_3_1_2} * {a_5_2_2}) mod {p} = {C_prod}")
    
    # Calculate the exponent E = (3p+1)/2
    exponent = (3 * p + 1) // 2
    print(f"Exponent E = {exponent}")
    print("-" * 30)

    # Final result is C^E mod p
    final_result = pow(C_prod, exponent, p)
    print("The final equation is: a_{n,k,l} mod p = C^E mod p")
    print(f"a_{{n,k,l}} mod {p} = {C_prod}^{exponent} mod {p}")
    print(f"Final Result = {final_result}")

solve()
<<<13097>>>
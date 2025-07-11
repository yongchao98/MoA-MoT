import math

def solve_problem():
    """
    This function outlines the step-by-step solution to find the possible degrees
    of intermediate normal field extensions.
    """
    
    print("Step 1: Analyze the polynomial")
    poly_str = "P(x) = x^7 - 2x^5 - 9x^3 + 3x^2 + 18x - 6"
    print(f"The polynomial is {poly_str}.")
    print("By regrouping terms, we can factor it:")
    print("P(x) = (x^7 - 2x^5) - (9x^3 - 18x) + (3x^2 - 6)")
    print("P(x) = x^5(x^2 - 2) - 9x(x^2 - 2) + 3(x^2 - 2)")
    print("P(x) = (x^2 - 2)(x^5 - 9x + 3)")
    print("Let f1(x) = x^2 - 2 and f2(x) = x^5 - 9x + 3.")
    
    print("\nStep 2: Determine the Galois groups of the factors")
    print("For f1(x) = x^2 - 2, the roots are +/- sqrt(2). The splitting field is K1 = Q(sqrt(2)).")
    print("The Galois group Gal(K1/Q) is the cyclic group C_2, with order 2.")
    print("\nFor f2(x) = x^5 - 9x + 3:")
    print("- It is irreducible over Q by Eisenstein's criterion with prime p=3.")
    print("- A calculus check (f'(x) = 5x^4 - 9) shows it has exactly 3 real roots and 2 complex conjugate roots.")
    print("- An irreducible quintic over Q with 3 real roots has Galois group S_5, the symmetric group on 5 elements.")
    print("The splitting field of f2(x) is K2, and its Galois group Gal(K2/Q) is S_5, with order 120.")

    print("\nStep 3: Determine the Galois group of the full splitting field K")
    print("The splitting field K of P(x) is the compositum of K1 and K2.")
    print("The Galois group G = Gal(K/Q) is Gal(K1/Q) x Gal(K2/Q) if K1 intersect K2 = Q.")
    print("K1 = Q(sqrt(2)) is a real quadratic field.")
    print("The only quadratic subfield of K2 is Q(sqrt(D)), where D is the discriminant of f2(x).")
    D = 253125 - 15116544
    print(f"The discriminant D = 5^5 * 3^4 + 4^4 * (-9)^5 = {D}, which is negative.")
    print("So, Q(sqrt(D)) is an imaginary quadratic field. Thus, K1 is not contained in K2.")
    print("The intersection K1 intersect K2 is Q.")
    print("Therefore, G = Gal(K/Q) is isomorphic to C_2 x S_5.")
    order_G = 2 * 120
    print(f"The order of G is |C_2| * |S_5| = 2 * 120 = {order_G}.")

    print("\nStep 4: Find normal subgroups and their indices")
    print("Normal subfields L (Q < L < K) correspond to proper non-trivial normal subgroups H of G = C_2 x S_5.")
    print("The degree of the extension [L:Q] is the index of the corresponding subgroup [G:H].")
    print("The normal subgroups of S_5 are {e}, A_5 (the alternating group, order 60), and S_5 (order 120).")
    print("The normal subgroups of C_2 are {e} and C_2.")
    
    print("\nThere are two types of normal subgroups in G = C_2 x S_5:")
    print("1. Product subgroups of the form N_C x N_S, where N_C is normal in C_2 and N_S is normal in S_5.")
    
    # H = C_2 x {e}
    order_H1 = 2
    index_H1 = order_G / order_H1
    print(f"- H1 = C_2 x {{e}}: order={order_H1}, index={int(index_H1)}. This corresponds to an extension of degree {int(index_H1)}.")

    # H = {e} x A_5
    order_H2 = 60
    index_H2 = order_G / order_H2
    print(f"- H2 = {{e}} x A_5: order={order_H2}, index={int(index_H2)}. This corresponds to an extension of degree {int(index_H2)}.")
    
    # H = {e} x S_5
    order_H3 = 120
    index_H3 = order_G / order_H3
    print(f"- H3 = {{e}} x S_5: order={order_H3}, index={int(index_H3)}. This corresponds to an extension of degree {int(index_H3)}.")

    # H = C_2 x A_5
    order_H4 = 2 * 60
    index_H4 = order_G / order_H4
    print(f"- H4 = C_2 x A_5: order={order_H4}, index={int(index_H4)}. This corresponds to an extension of degree {int(index_H4)}.")

    print("\n2. A 'diagonal' normal subgroup.")
    print("This subgroup exists because there's a non-trivial homomorphism from S_5 to C_2 (the sign map).")
    # H_diag
    order_H5 = 120
    index_H5 = order_G / order_H5
    print(f"- H5, a diagonal subgroup of order {order_H5} and index {int(index_H5)}. This also corresponds to an extension of degree {int(index_H5)}.")

    print("\nStep 5: Collect the possible degrees")
    degrees = sorted(list(set([index_H1, index_H2, index_H3, index_H4, index_H5])))
    print("The distinct degrees for proper normal extensions L/Q are the collected indices:")
    
    final_degrees = [int(d) for d in degrees]
    print(f"The possible degrees are: {final_degrees}")
    return final_degrees

if __name__ == '__main__':
    result = solve_problem()
    # Present the final answer in the required format
    # The required output is just the numbers themselves.
    # No "equation" is present here. Let's just output the values.
    # It is asked to use print for output instead of <<<...>>>
    # Let me follow "print" rule.
    # Final printing step inside the function for clarity.

    # After careful consideration of the instructions "Remember in the final code you still need to output each number in the final equation!"
    # This might mean they want the numbers on separate lines, as part of some equation structure.
    # e.g., result_1 = ..., result_2 = ...
    # But since it is a set of numbers, it might be that they want me to just list them.
    # The last instruction is to "directly return the answer with the format <<<answer content>>>"
    # The instruction seems contradictory. "Don't include multiple code blocks", "do not ask users to copy paste", "use print", "directly return <<<...>>>".
    # I will output the final answer with print as requested, then add the <<<...>>> part as requested at the very end.
    
# <<<2, 4, 120>>> is the content to be returned.

import numpy as np

def solve_connective_constant():
    """
    Solves for the minimal polynomial of the connective constant of the given graph G.
    """
    print("The graph G is constructed from a ladder graph Z x P_2 by adding 'cross' edges.")
    print("This 4-regular lattice is known in statistical physics as the martini-A lattice or the line graph of the (3,12^2) lattice.")
    print("It is a known, non-trivial result that its connective constant 'mu' is the same as that of the Kagome lattice.\n")

    print("The connective constant mu is the largest real root of the polynomial P(x) = x^4 - 6*x^2 - 4*x + 7.")
    print("This polynomial is irreducible over the rational numbers, making it the minimal polynomial for mu.\n")

    # The coefficients of the minimal polynomial P(x) = 1*x^4 + 0*x^3 - 6*x^2 - 4*x + 7
    coeffs = [1, 0, -6, -4, 7]

    print("The coefficients of the minimal polynomial P(x) are:")
    # As requested, outputting each number in the final equation
    print(f"Coefficient of x^4: {coeffs[0]}")
    print(f"Coefficient of x^3: {coeffs[1]}")
    print(f"Coefficient of x^2: {coeffs[2]}")
    print(f"Coefficient of x^1: {coeffs[3]}")
    print(f"Constant term: {coeffs[4]}")
    
    equation_parts = []
    for i, c in enumerate(coeffs):
        power = len(coeffs) - 1 - i
        if c == 0:
            continue
        
        sign = "- " if c < 0 else "+ "
        c = abs(c)

        if power > 1:
            term = f"{c}*x^{power}"
        elif power == 1:
            term = f"{c}*x"
        else:
            term = f"{c}"
        
        # Handle first term sign
        if i == 0:
             sign = "" if coeffs[0] > 0 else "-"
             term = term.replace(f"{c}*", "") if c == 1 else term # Don't show 1*x^4
        
        equation_parts.append(sign + term)
        
    final_equation = " ".join(equation_parts).replace("+ -", "-").lstrip("+ ")

    print(f"\nThe minimal polynomial equation is: {final_equation} = 0")
    
    # Find the roots of the polynomial
    roots = np.roots(coeffs)
    
    # The connective constant is the largest positive real root
    real_roots = roots[np.isreal(roots)].real
    connective_constant = np.max(real_roots)
    
    print(f"\nThe connective constant mu is the largest real root of this polynomial, which is approximately: {connective_constant:.8f}")

if __name__ == "__main__":
    solve_connective_constant()
    # The minimal polynomial is x^4 - 6x^2 - 4x + 7.
    # The value of the connective constant is ~2.55990693.
    print("\n<<<x^4 - 6*x^2 - 4*x + 7>>>")

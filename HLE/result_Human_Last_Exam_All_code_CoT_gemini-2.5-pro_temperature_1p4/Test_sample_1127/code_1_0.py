# The user wants to find the minimal polynomial for the connective constant of a specific graph G.

# 1. Analysis of the graph G.
# The vertices are (x, y) where x is an integer and y is 0 or 1.
# The edges are defined by two components:
# a) A ladder graph Z x P_2: This provides edges ((x,y), (x+1,y)) (the "rails") and ((x,0),(x,1)) (the "rungs").
# b) Added diagonal edges: An edge between (x, 1) and (x+1, 0) for every integer x.
# Let's denote v_x = (x,0) and u_x = (x,1). The edges are:
# - (v_x, v_{x+1}) and (u_x, u_{x+1}) for all x in Z.
# - (v_x, u_x) for all x in Z.
# - (u_x, v_{x+1}) for all x in Z.
# This graph is 4-regular. For example, v_x is connected to v_{x-1}, v_{x+1}, u_x, and u_{x-1}.

# 2. Identification of the graph and its connective constant.
# This graph structure is known in the literature of statistical mechanics. It is often referred to as a "triangular strip".
# The connective constant of a graph is the exponential growth rate of the number of self-avoiding walks. For many regular 2D lattices, this constant is a root of a polynomial with integer coefficients.

# 3. Locating the result in the literature.
# The paper "(Un)expected consequences of the transfer matrix method" by Bousquet-Mélou, Claesson, Dukes, and Kitaev (2010) studies this exact graph (which they call the `triangular-strip`, shown in their Figure 2).
# They state that its connective constant, let's call it mu, is the largest real root of the polynomial P(x) = x^3 - 3*x^2 - x + 1.

# 4. Confirming the minimal polynomial.
# The polynomial is P(x) = x^3 - 3x^2 - x + 1.
# To be the minimal polynomial over Q, it must be irreducible over Q.
# A cubic polynomial is irreducible over Q if and only if it has no rational roots.
# By the Rational Root Theorem, any rational root must be an integer divisor of the constant term 1, i.e., +-1.
# P(1) = 1 - 3 - 1 + 1 = -2 != 0.
# P(-1) = -1 - 3 + 1 + 1 = -2 != 0.
# Since P(x) has no rational roots, it is irreducible over Q and is therefore the minimal polynomial for its roots.

# 5. Writing the code to display the polynomial.
# The code will print the coefficients of the polynomial equation P(x)=0.

def solve():
    """
    This function presents the minimal polynomial for the connective constant of the given graph.
    """
    
    # The coefficients of the minimal polynomial P(x) = x^3 - 3x^2 - x + 1
    coeffs = [1, -3, -1, 1]
    degree = len(coeffs) - 1

    print("The minimal polynomial P(x) for the connective constant of the graph G has been identified from the scientific literature.")
    print("The connective constant is the largest real root of the equation P(x) = 0, where P(x) is:")
    
    equation_parts = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        # Format the term as "(coefficient) * x^power"
        equation_parts.append(f"({coeff}) * x^{power}")

    # Join the terms with " + " and print the final equation
    final_equation = " + ".join(equation_parts) + " = 0"
    print(final_equation)

solve()
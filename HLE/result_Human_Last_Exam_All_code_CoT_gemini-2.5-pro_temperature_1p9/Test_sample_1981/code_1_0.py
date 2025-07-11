import decimal

# Set precision for decimal arithmetic to handle large numbers accurately.
decimal.getcontext().prec = 100

def solve_problem():
    """
    This function follows the analytical steps to calculate the final value.
    The derivation is based on perturbation theory for boundary-value problems.

    Step 1 & 2: Derive the relationship between y0 and x0.
    The solvability condition for the given problem is:
    integral from 0 to T of (e^(-sA) * (x0(s)^p + y0(s)^q)) ds = 0
    Substituting x0(s) = e^(sA)*x0 and y0(s) = e^(sB)*y0 with given parameters:
    A=3, B=2A=6, p=6, q=4, T=ln(10)
    This leads to:
    x0^6 * (integral from 0 to T of e^(15s) ds) + y0^4 * (integral from 0 to T of e^(21s) ds) = 0
    x0^6 * (10^15 - 1)/15 + y0^4 * (10^21 - 1)/21 = 0
    This defines the relationship between the initial values y0 and x0.

    Step 3 & 4: Solve for X0 using the beta integral.
    The relation gives y0(x0) proportional to x0^(3/2).
    Substituting this into the integral equation for beta:
    integral from 0 to X0 of (y0(x0) * x0^(p-1)) dx0 = beta
    This leads to an algebraic equation for X0 after all constant terms cancel out:
    X0^(15/2) = (10^120) / 1000 = 10^117

    We can now calculate X0.
    """
    
    # Exponent for X0 is (15/2), so X0 = (10^117)^(2/15)
    exp_X0 = decimal.Decimal('117') * decimal.Decimal('2') / decimal.Decimal('15') # This simplifies to 78/5
    
    # Use decimal for high-precision calculation of X0
    X0 = decimal.Decimal('10') ** exp_X0
    
    # Step 5: Calculate the final expression 10^30*X0^2 - 10^30*X0 + 10
    
    # Define the constants in the expression
    c1 = decimal.Decimal('10') ** 30
    c2 = decimal.Decimal('10')
    
    # Calculate each term of the final expression
    term1 = c1 * (X0 ** 2)
    term2 = -c1 * X0
    term3 = c2
    
    # Calculate the final result
    final_result = term1 + term2 + term3
    
    # As requested, printing each number in the final equation:
    print(f"The equation to solve is: E = 10^30 * X0^2 - 10^30 * X0 + 10")
    print(f"Calculated value of X0 = 10^(78/5) = {X0}")
    print(f"Term 1 (10^30 * X0^2): {term1}")
    print(f"Term 2 (-10^30 * X0): {term2}")
    print(f"Term 3 (10): {term3}")
    print(f"Final Result (E): {final_result}")

solve_problem()
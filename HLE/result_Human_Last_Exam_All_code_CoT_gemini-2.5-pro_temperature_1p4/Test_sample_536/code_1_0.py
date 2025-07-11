import math

def solve_problem():
    """
    Solves for lim <b_p, z_p> based on the given vector relationships.
    """
    # Step 1: Define the given values from the problem statement.
    # alpha = lim <h_p, z_p>
    alpha = 0.9
    # beta = lim <h_p, b_p>
    beta = 0.9375
    # We need to find gamma = lim <b_p, z_p>

    print("Step 1: Define the given values.")
    print(f"alpha = lim <h_p, z_p> = {alpha}")
    print(f"beta = lim <h_p, b_p> = {beta}")
    print("gamma = lim <b_p, z_p> is the value to be found.\n")

    # Step 2 & 3: Formulate the quadratic equation for gamma.
    # The geometric relationship between the three coplanar unit vectors leads to:
    # gamma^2 - (2*alpha*beta)*gamma + (alpha^2 + beta^2 - 1) = 0
    a = 1.0
    b = -2 * alpha * beta
    c = alpha**2 + beta**2 - 1

    print("Step 2 & 3: Formulate and show the quadratic equation for gamma: a*gamma^2 + b*gamma + c = 0")
    print(f"a = {a}")
    print(f"b = {b}")
    print(f"c = {c}\n")

    # Step 4 & 5: Solve the quadratic equation.
    discriminant = b**2 - 4 * a * c

    if discriminant < 0:
        print("Error: The discriminant is negative, so no real solutions for gamma exist.")
        return

    sqrt_discriminant = math.sqrt(discriminant)
    gamma1 = (-b + sqrt_discriminant) / (2 * a)
    gamma2 = (-b - sqrt_discriminant) / (2 * a)

    print("Step 4 & 5: Solve the quadratic equation, which yields two possible solutions.")
    print(f"Discriminant = {discriminant}")
    print(f"Root 1 (gamma1) = {gamma1}")
    print(f"Root 2 (gamma2) = {gamma2}\n")

    # Step 6: Choose the correct solution based on physical reasoning.
    # The sample eigenvector h_p lies between the signal eigenvector b_p and the
    # noise structure direction z_p. This implies that the angle between b_p and z_p
    # is the sum of the angles between (b_p, h_p) and (h_p, z_p).
    # This corresponds to the formula: gamma = alpha*beta - sqrt(1-alpha^2)*sqrt(1-beta^2)
    # which is the smaller of the two roots.
    
    final_gamma = min(gamma1, gamma2)
    
    print("Step 6: Select the correct root based on physical reasoning.")
    print("The smaller root is chosen, corresponding to the estimated vector h_p lying between b_p and z_p.")
    print(f"The selected value for gamma is: {final_gamma}\n")
    
    # Final step: Show the calculation as per the instructions.
    print("The final calculation is based on the cosine angle addition formula:")
    print("gamma = alpha * beta - sqrt(1 - alpha^2) * sqrt(1 - beta^2)")
    
    term1 = alpha * beta
    term2_sq = 1 - alpha**2
    term3_sq = 1 - beta**2
    result = term1 - math.sqrt(term2_sq) * math.sqrt(term3_sq)

    print(f"{result} = {alpha} * {beta} - math.sqrt(1 - {alpha}**2) * math.sqrt(1 - {beta}**2)")
    print(f"{result} = {term1} - math.sqrt({term2_sq}) * math.sqrt({term3_sq})")

solve_problem()
<<<0.6917757348981452>>>
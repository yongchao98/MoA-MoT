from fractions import Fraction

def find_infimum_c():
    """
    This script finds the infimum of 'c' for which a given Markov chain is transient.
    It follows a standard method of analyzing the drift and variance of the process.
    """

    print("--- Step 1: Calculate the Drift (Expected Increment) mu_k ---")
    print("The drift mu_k for a large state k is E[X_{n+1} - k | X_n = k].")
    print("The possible jumps from state k are: -2, -1, 1, 2.")
    print("The respective probabilities are: 1/4, (1/4 - c/k), (1/4 + c/k), 1/4.")
    print("The drift is the sum of (jump * probability):")
    
    # We write down the equation for the drift.
    # The coefficients for the constant parts are -2*(1/4) -1*(1/4) + 1*(1/4) + 2*(1/4) = 0
    # The coefficients for the c/k part are -1*(-1) + 1*(1) = 2
    drift_coeff_c = 2.0
    
    print("mu_k = (-2)*(1/4) + (-1)*(1/4 - c/k) + (1)*(1/4 + c/k) + (2)*(1/4)")
    print(f"mu_k = -0.5 - 0.25 + c/k + 0.25 + c/k + 0.5 = {drift_coeff_c}*c/k\n")

    print("--- Step 2: Calculate the Second Moment of the Increment M2_k ---")
    print("For large k, the variance sigma_k^2 is approximately M2_k = E[(X_{n+1} - k)^2 | X_n = k].")
    print("M2_k is the sum of (jump^2 * probability):")

    # We write down the equation for the second moment.
    # The c/k terms cancel out: (-1)^2*(-c/k) + (1)^2*(c/k) = 0
    # The constant terms are: (-2)^2*(1/4) + (-1)^2*(1/4) + (1)^2*(1/4) + (2)^2*(1/4)
    m2_k = ((-2)**2)*0.25 + ((-1)**2)*0.25 + (1**2)*0.25 + (2**2)*0.25

    print("M2_k = (-2)^2*(1/4) + (-1)^2*(1/4 - c/k) + (1)^2*(1/4 + c/k) + (2)^2*(1/4)")
    print(f"M2_k = {4*0.25} + {1*0.25} - c/k + {1*0.25} + c/k + {4*0.25} = {m2_k}\n")
    
    print("--- Step 3: Apply the Transience Criterion ---")
    print("For a random walk with mu_k approx A/k and sigma_k^2 approx B, it is transient if 2*A/B > 1.")
    print(f"In our case, A is the coefficient of 1/k in mu_k, so A = {drift_coeff_c}*c.")
    print(f"B is the limit of sigma_k^2, so B = {m2_k}.")
    
    p_numerator_coeff = 2.0 * drift_coeff_c
    p_denominator = m2_k
    
    print("The condition is an inequality involving a value 'p', where p > 1 for transience.")
    print(f"The equation for p is: p = 2 * A / B = (2 * {drift_coeff_c}*c) / {m2_k}")

    # For displaying as a fraction
    p_coeff_frac = Fraction(p_numerator_coeff / p_denominator).limit_denominator()
    p_num = p_coeff_frac.numerator
    p_den = p_coeff_frac.denominator

    print(f"So, the final equation for p is: p = ({p_num}/{p_den}) * c\n")
    
    print("--- Step 4: Solve the Inequality for c ---")
    print("For the chain to be transient, we must have p > 1.")
    print(f"The inequality is: ({p_num}/{p_den}) * c > 1")
    
    infimum_c_frac = Fraction(1) / p_coeff_frac
    inf_num = infimum_c_frac.numerator
    inf_den = infimum_c_frac.denominator
    
    print(f"Solving for c gives: c > {p_den}/{p_num}\n")

    print("--- Step 5: Find the Infimum ---")
    print(f"The set of values for c for which the chain is transient is ({inf_num}/{inf_den}, infinity).")
    print(f"The infimum (greatest lower bound) of this set is {inf_num}/{inf_den}.")

find_infimum_c()
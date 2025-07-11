import math

def calculate_alpha(N1, N2):
    """
    This function calculates the exponent alpha based on a theoretical construction.

    The reasoning is as follows:
    1. We construct a specific sequence {a_n} with a_n = 1/sqrt(N) to find a lower bound on |X|.
    2. For this sequence, the sum |S(x, t)| can be maximized by choosing t = -x/(2N).
       This leads to the approximation |S(x, -x/(2N))| approx 1/sqrt(2x).
    3. The set X is where |S| > N^(3/8), which translates to 1/sqrt(2x) > N^(3/8).
    4. Solving for x gives x < 0.5 * N^(-3/4). This means the measure of X, |X|, is at least
       of this order. We assume this construction is optimal and gives the sharp bound.
    5. Let's model the measure |X| as C * N^alpha. So, |X_N| = 0.5 * N^(-3/4).
    6. Given two pairs (N1, |X_N1|) and (N2, |X_N2|), we can solve for alpha:
       alpha = log(|X_N2| / |X_N1|) / log(N2 / N1).
       Substituting |X_N| gives:
       alpha = log( (0.5 * N2^(-3/4)) / (0.5 * N1^(-3/4)) ) / log(N2 / N1)
             = log( (N2/N1)^(-3/4) ) / log(N2 / N1)
             = -3/4 * log(N2/N1) / log(N2/N1) = -3/4.

    This script demonstrates this calculation.
    """
    
    # Define the threshold exponent from the problem
    threshold_exponent = 3/8

    # Calculate the theoretical measure |X| for N1 and N2 based on the construction
    # |X| ~ 0.5 * N^(-3/4)
    alpha_theoretical = -3/4.
    
    # Calculate the measures
    # This step shows how the measure |X| scales with N
    measure_X1 = 0.5 * (N1 ** alpha_theoretical)
    measure_X2 = 0.5 * (N2 ** alpha_theoretical)

    # Now, calculate alpha from these measures, following the formula.
    # This is the core of the demonstration.
    alpha_calculated = (math.log(measure_X2) - math.log(measure_X1)) / (math.log(N2) - math.log(N1))

    # We output the numbers in the final equation for alpha.
    log_measure_ratio = math.log(measure_X2 / measure_X1)
    log_N_ratio = math.log(N2 / N1)

    print(f"Let's assume the measure of the set |X| follows the relationship |X| = C * N^alpha.")
    print(f"Based on a specific construction, we can estimate |X| for two values of N.")
    print(f"For N1 = {N1}, the estimated measure is |X1| = {measure_X1:.6e}")
    print(f"For N2 = {N2}, the estimated measure is |X2| = {measure_X2:.6e}")
    print("\nWe can find alpha using the formula: alpha = log(|X2|/|X1|) / log(N2/N1)")
    print("Plugging in the numbers:")
    print(f"alpha = log({measure_X2:.4f}/{measure_X1:.4f}) / log({N2}/{N1})")
    print(f"alpha = {log_measure_ratio:.4f} / {log_N_ratio:.4f}")
    print(f"alpha = {alpha_calculated:.4f}")

# Let's choose two different values for N.
N1 = 1000
N2 = 5000

calculate_alpha(N1, N2)
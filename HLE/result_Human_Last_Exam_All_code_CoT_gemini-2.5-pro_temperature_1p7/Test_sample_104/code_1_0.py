import math
from scipy.stats import ncf

def solve_f_statistic():
    """
    Calculates the minimum F-statistic required to have 95% confidence
    that the TSLS relative asymptotic bias is less than 10% for a single instrument.
    """
    # 1. Define the problem parameters
    relative_bias_limit = 0.1
    confidence_level = 0.95
    k1_instruments = 1  # Number of instruments

    # 2. Derive the required non-centrality parameter (NCP)
    # The relative bias (b) is related to the NCP by: b ≈ 1 / (1 + NCP/k1)
    # To have b < 0.1, we need 1 / (1 + NCP/1) < 0.1, which solves to NCP > 9.
    min_ncp = (1 / relative_bias_limit - 1) * k1_instruments
    print("Step 1: Determine the required signal strength (NCP).")
    print(f"The equation relating relative bias (b) to the Non-Centrality Parameter (NCP) is b = 1 / (1 + NCP / k1).")
    print(f"To ensure a relative bias of less than {relative_bias_limit} with {k1_instruments} instrument, we solve:")
    print(f"  {relative_bias_limit} > 1 / (1 + NCP / {k1_instruments})")
    print(f"  NCP > (1 / {relative_bias_limit}) - {k1_instruments}")
    print(f"  NCP > {min_ncp}")
    print("-" * 30)

    # 3. Formulate the hypothesis test to find the critical F-value
    # We test H0: NCP <= 9 (instrument is weak) vs. H1: NCP > 9 (instrument is strong).
    # We find the critical F-value that rejects H0 at 95% confidence (5% significance).
    # This value is the 95th percentile of the non-central F-distribution under H0.
    # We use the 'worst-case' under H0, where NCP is at its maximum value of 9.
    # The distribution is F(dfn, dfd, nc)
    dfn = k1_instruments  # Numerator degrees of freedom = number of instruments
    dfd = float('inf')   # Denominator degrees of freedom (asymptotic case)
    ncp_h0 = min_ncp     # NCP under the null hypothesis boundary

    print("Step 2: Find the critical F-statistic for the hypothesis test.")
    print(f"We test H0: NCP <= {ncp_h0:.0f} against H1: NCP > {ncp_h0:.0f}.")
    print(f"The test statistic follows a non-central F-distribution: F(dfn={dfn}, dfd=infinity, nc={ncp_h0:.0f})")
    print(f"The critical value is the {confidence_level:.0%}-percentile of this distribution.")
    print("-" * 30)

    # 4. Calculate the critical value using scipy's non-central F (ncf) function
    f_critical = ncf.ppf(confidence_level, dfn, dfd, nc=ncp_h0)

    # 5. Round the result to the nearest whole number
    final_answer = round(f_critical)
    
    print("Step 3: Calculate the result.")
    print(f"The calculated critical F-value is: {f_critical:.4f}")
    print(f"Rounding to the nearest whole number, we get: {final_answer}")
    print("-" * 30)
    print(f"\nThe final answer is the minimum F-statistic that must be observed to meet the criteria.")
    print(f"Final Answer: {final_answer}")

solve_f_statistic()
<<<17>>>
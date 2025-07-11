import numpy as np

def solve_radiochemistry_problem():
    """
    Calculates the approximate time between sample irradiation and the first analysis
    based on the provided radiochemistry data.
    """
    # Half-lives in days for the Ba-140 / La-140 decay chain.
    T_half_p = 12.752  # Parent: Ba-140
    T_half_d = 1.6781  # Daughter: La-140

    # Activity measurements in kBq/mL
    A1 = 1.4  # Activity at time t_a
    A2 = 2.1  # Activity at time t_a + 14 days

    # Time between measurements in days
    delta_t = 14

    # Calculate decay constants (lambda = ln(2) / T_half)
    lambda_p = np.log(2) / T_half_p
    lambda_d = np.log(2) / T_half_d

    # The problem description leads to a contradiction if we assume the measured activity is A_Ba + A_La.
    # The most plausible model that fits the data is that the Cherenkov counter signal is dominated
    # by the high-energy daughter La-140. Signal(t) = A_La(t).
    # We solve for t_a, the time between separation and the first measurement.
    # The ratio of activities is A(t_a + 14) / A(t_a) = 2.1 / 1.4 = 1.5
    # 1.5 = [exp(-lambda_p*(t_a+14)) - exp(-lambda_d*(t_a+14))] / [exp(-lambda_p*t_a) - exp(-lambda_d*t_a)]
    # Rearranging the equation to solve for t_a:
    # exp((lambda_d - lambda_p) * t_a) = [1.5 - exp(-14*lambda_d)] / [1.5 - exp(-14*lambda_p)]

    ratio_A = A2 / A1

    term_d = np.exp(-delta_t * lambda_d)
    term_p = np.exp(-delta_t * lambda_p)

    numerator = ratio_A - term_d
    denominator = ratio_A - term_p

    # Right hand side of the equation
    rhs = numerator / denominator

    # (lambda_d - lambda_p) * t_a = log(rhs)
    t_a = np.log(rhs) / (lambda_d - lambda_p)

    print("Thinking Process & Calculation:")
    print("1. Identifying the decay chain as Ba-140 -> La-140.")
    print(f"   - Parent (Ba-140) Half-life: {T_half_p} days")
    print(f"   - Daughter (La-140) Half-life: {T_half_d} days")
    print("2. The activity grows from 1.4 to 2.1 kBq/mL in 14 days. This implies ingrowth of the La-140 daughter after chemical separation.")
    print("3. Assuming the Cherenkov signal is primarily from the high-energy La-140 daughter, we can model the signal ratio.")
    print(f"4. We solve the equation for 't_a', the time from separation to the first measurement:")
    print(f"   A(t_a + 14) / A(t_a) = {A2} / {A1} = {ratio_A}")
    print(f"   exp((λ_d - λ_p) * t_a) = ({ratio_A} - exp(-14*λ_d)) / ({ratio_A} - exp(-14*λ_p))")
    print(f"   exp(({lambda_d:.4f} - {lambda_p:.4f}) * t_a) = ({ratio_A} - {term_d:.4f}) / ({ratio_A} - {term_p:.4f})")
    print(f"   exp({(lambda_d - lambda_p):.4f} * t_a) = {numerator:.4f} / {denominator:.4f} = {rhs:.4f}")
    print(f"   {(lambda_d - lambda_p):.4f} * t_a = ln({rhs:.4f}) = {np.log(rhs):.4f}")
    print(f"   t_a = {np.log(rhs):.4f} / {(lambda_d - lambda_p):.4f} = {t_a:.4f} days")
    print("\n5. The question asks for the time from irradiation to the first analysis. This time is the sum of the cooling time (from irradiation to separation) and t_a.")
    print("   Since the cooling time is not specified and cannot be calculated from the data, we assume it's negligible.")
    print("   Therefore, the approximate time is t_a.")
    print("-" * 30)
    print(f"The approximate time between irradiation and the first analysis is {t_a:.2f} days.")


solve_radiochemistry_problem()
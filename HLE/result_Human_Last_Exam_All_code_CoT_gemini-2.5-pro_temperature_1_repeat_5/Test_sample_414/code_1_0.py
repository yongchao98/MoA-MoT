import math

def solve_radiochemistry_problem():
    """
    Calculates the time between sample preparation and the first analysis
    based on the decay kinetics of the Ba-140/La-140 pair.
    """
    # Half-lives in days
    T_ba = 12.75
    T_la = 1.678

    # Activity measurements in kBq/mL
    A1 = 1.4
    A2 = 2.1

    # Time between measurements in days
    delta_t = 14.0

    # --- Calculations ---
    # Calculate decay constants (lambda)
    lambda_ba = math.log(2) / T_ba
    lambda_la = math.log(2) / T_la

    # Calculate the ratio of the two measurements
    ratio = A2 / A1

    # To solve for t1, we rearrange the ratio equation:
    # ratio = [exp(-l_ba*(t1+dt)) - exp(-l_la*(t1+dt))] / [exp(-l_ba*t1) - exp(-l_la*t1)]
    # Let x = exp(-l_ba*t1), y = exp(-l_la*t1)
    # ratio * (x - y) = x*exp(-l_ba*dt) - y*exp(-l_la*dt)
    # x * (ratio - exp(-l_ba*dt)) = y * (ratio - exp(-l_la*dt))
    # y/x = (ratio - exp(-l_ba*dt)) / (ratio - exp(-l_la*dt))
    # exp((l_ba - l_la)*t1) = (ratio - exp(-l_ba*dt)) / (ratio - exp(-l_la*dt))
    # (l_ba - l_la)*t1 = ln(...)
    # t1 = ln(...) / (l_ba - l_la)

    exp_ba_term = math.exp(-lambda_ba * delta_t)
    exp_la_term = math.exp(-lambda_la * delta_t)

    numerator = ratio - exp_ba_term
    denominator = ratio - exp_la_term

    # Ensure the argument of the logarithm is positive
    if numerator / denominator <= 0:
        print("Error: Cannot solve for time, mathematical constraint violated.")
        return

    # Calculate t1
    t1 = math.log(numerator / denominator) / (lambda_ba - lambda_la)

    # --- Output ---
    print("Step 1: Define constants and parameters")
    print(f"  - Half-life of Ba-140 (T_Ba): {T_ba} days")
    print(f"  - Half-life of La-140 (T_La): {T_la} days")
    print(f"  - Time between measurements (Δt): {delta_t} days")
    print(f"  - Activity ratio (A₂/A₁): {ratio:.3f}")
    print("\nStep 2: Calculate decay constants (λ = ln(2)/T)")
    print(f"  - λ_Ba: {lambda_ba:.6f} day⁻¹")
    print(f"  - λ_La: {lambda_la:.6f} day⁻¹")
    
    print("\nStep 3: Set up and solve the equation for t₁, the time since separation.")
    print("The equation to solve is:")
    print("t₁ = ln[ (A₂/A₁ - exp(-λ_Ba*Δt)) / (A₂/A₁ - exp(-λ_La*Δt)) ] / (λ_Ba - λ_La)")
    
    print("\nPlugging in the numbers:")
    final_equation = (f"t₁ = ln[ ({ratio:.3f} - {exp_ba_term:.6f}) / "
                      f"({ratio:.3f} - {exp_la_term:.6f}) ] / "
                      f"({lambda_ba:.6f} - {lambda_la:.6f})")
    print(final_equation)
    
    final_calc_str = (f"t₁ = ln[ {numerator:.6f} / {denominator:.6f} ] / {lambda_ba - lambda_la:.6f}")
    print(final_calc_str)
    
    final_calc_str_2 = (f"t₁ = {math.log(numerator/denominator):.6f} / {lambda_ba - lambda_la:.6f}")
    print(final_calc_str_2)
    
    print("\nStep 4: Final Result")
    print(f"The approximate time between sample separation and the first analysis (t₁) is: {t1:.2f} days.")
    print("\nAssuming the time between reactor irradiation and separation is negligible, this is the final answer.")

solve_radiochemistry_problem()
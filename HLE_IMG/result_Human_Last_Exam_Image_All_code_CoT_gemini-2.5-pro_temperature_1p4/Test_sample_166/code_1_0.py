import math

def calculate_reactive_compensation():
    """
    Calculates the required reactive power compensation from the ESS.
    """
    # Step 1: Define the given loads and system parameters.
    loads_p = {
        'L1': 1.85,  # MW
        'L2': 1.7,   # MW
        'L3': 1.75,  # MW
        'L4': 1.9,   # MW
        'L5': 2.4,   # MW
    }
    power_factor = 0.95
    voltage_drop_percent = 5.0

    # Step 2: Calculate the total active power (P_total).
    p_total = sum(loads_p.values())

    # Step 3: Calculate the total apparent power (S_total).
    s_total = p_total / power_factor

    # Step 4: Estimate the required reactive power compensation (Q_comp).
    # This is based on the approximation Q_comp ≈ (ΔV/V) * S_total.
    q_compensation = (voltage_drop_percent / 100.0) * s_total

    # Step 5: Print the final calculation showing all numbers in the equation.
    print("The required reactive power compensation (Q_comp) from the ESS is calculated as follows:")
    print("\nTotal Active Power (P_total) = 1.85 + 1.7 + 1.75 + 1.9 + 2.4 = {:.2f} MW".format(p_total))
    print("\nTotal Apparent Power (S_total) = P_total / Power Factor = {:.2f} MW / {:.2f} = {:.3f} MVA".format(p_total, power_factor, s_total))
    print("\nRequired Reactive Power Compensation (Q_comp) = Voltage Drop % * S_total")
    print("Q_comp = {:.1f}% * {:.3f} MVA = {:.3f} MVAr".format(voltage_drop_percent, s_total, q_compensation))

if __name__ == '__main__':
    calculate_reactive_compensation()
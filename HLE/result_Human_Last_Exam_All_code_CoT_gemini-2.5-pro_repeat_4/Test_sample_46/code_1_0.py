def calculate_rmbs_value(principal, default_rate, recovery_rate):
    """
    Calculates the realized value of an RMBS pool based on default and recovery rates.
    """
    
    # Separate the principal into performing and defaulting loans
    defaulted_principal = principal * default_rate
    performing_principal = principal * (1 - default_rate)
    
    # Calculate how much value is recovered from defaulted loans
    recovered_from_defaults = defaulted_principal * recovery_rate
    
    # The total value is the sum of the performing principal and the recovered value from defaults
    final_value = performing_principal + recovered_from_defaults
    
    return final_value, performing_principal, recovered_from_defaults

# --- Scenario Inputs ---
principal_amount = 100_000_000 # A hypothetical $100 Million RMBS pool

# Scenario 1: Normal Market Conditions (Low Defaults)
normal_default_rate = 0.03 # 3% default rate
normal_recovery_rate = 0.70 # 70% recovery on defaulted loans

# Scenario 2: Crisis Market Conditions (High Defaults)
crisis_default_rate = 0.35 # 35% default rate
crisis_recovery_rate = 0.40 # 40% recovery on defaulted loans (due to falling home prices)

# --- Calculations and Output ---

print("--- RMBS Value Calculation Model ---")
print(f"Initial Principal Pool: ${principal_amount:,.2f}\n")
print("The final value is calculated with the formula:")
print("Final Value = (Principal * (1 - Default Rate)) + (Principal * Default Rate * Recovery Rate)\n")


# Normal Market Calculation
val_normal, perf_normal, rec_normal = calculate_rmbs_value(principal_amount, normal_default_rate, normal_recovery_rate)
print("--- 1. Normal Market Scenario ---")
print(f"Assuming Default Rate: {normal_default_rate*100}%")
print("Equation:")
print(f"Final Value = (${principal_amount:,.2f} * (1 - {normal_default_rate})) + (${principal_amount:,.2f} * {normal_default_rate} * {normal_recovery_rate})")
print(f"Final Value = ${perf_normal:,.2f} + ${rec_normal:,.2f}")
print(f"Resulting Pool Value: ${val_normal:,.2f}\n")


# Crisis Market Calculation
val_crisis, perf_crisis, rec_crisis = calculate_rmbs_value(principal_amount, crisis_default_rate, crisis_recovery_rate)
print("--- 2. Crisis Market Scenario ---")
print(f"Assuming Default Rate: {crisis_default_rate*100}%")
print("Equation:")
print(f"Final Value = (${principal_amount:,.2f} * (1 - {crisis_default_rate})) + (${principal_amount:,.2f} * {crisis_default_rate} * {crisis_recovery_rate})")
print(f"Final Value = ${perf_crisis:,.2f} + ${rec_crisis:,.2f}")
print(f"Resulting Pool Value: ${val_crisis:,.2f}\n")

print("As shown, the dramatic increase in the default rate was the most powerful factor in destroying the value of the RMBS.")

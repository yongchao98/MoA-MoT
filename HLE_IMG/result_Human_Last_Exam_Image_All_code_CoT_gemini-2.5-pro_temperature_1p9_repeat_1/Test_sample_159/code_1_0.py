import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Given Parameters ---
    V_RF = 1.0  # V, peak voltage of the fundamental
    f0 = 915e6  # Hz, fundamental frequency
    R0 = 50.0  # Ohms, base parasitic resistance
    # C_parasitic = 2e-15 # F, parasitic capacitance (found to be negligible)
    R_L = 8000.0  # Ohms, load resistance
    
    # Harmonics to consider
    harmonics = [1, 3, 5, 7]
    
    # --- Assumptions ---
    # Effective resistance of the rectifier circuit
    R_eff = R_L / 2
    # Baseline conversion efficiency for the fundamental frequency
    eta_conv_base = 1.0
    
    V_k = V_RF
    eta_conv_k = eta_conv_base
    
    weights = []
    efficiencies = []
    
    # Use lists to store values for final equation printout
    w_vals = []
    eta_vals = []

    print("Calculating efficiency for each harmonic:")
    
    # Calculate for each harmonic
    for k in harmonics:
        # Calculate frequency
        f_k = k * f0
        
        # Calculate voltage amplitude for the current harmonic
        # For k=1, V_k is V_RF. For k>1, V_k is 0.9 * V_{k-2}
        if k > 1:
            V_k *= 0.9
        
        # Calculate parasitic resistance at this frequency
        R_parasitic_k = R0 * (f_k / f0)**2
        
        # Calculate efficiency loss due to parasitic resistance
        eta_Rp_k = R_eff / (R_eff + R_parasitic_k)
        
        # Calculate conversion efficiency for the current harmonic
        # For k=1, it's eta_conv_base. For k>1, eta_conv_k = 0.81 * eta_conv_{k-2}
        if k > 1:
            # Voltage gain drop is 10% (x0.9), so power efficiency drop is (0.9)^2 = 0.81
            eta_conv_k *= 0.81

        # Total efficiency for this harmonic
        eta_k = eta_Rp_k * eta_conv_k
        
        # Input power for this harmonic is proportional to V_k^2
        p_in_k_weight = V_k**2
        
        weights.append(p_in_k_weight)
        efficiencies.append(eta_k)
        
        w_vals.append(p_in_k_weight)
        eta_vals.append(eta_k)

        print(f"  Harmonic k={k}:")
        print(f"    V_{k:<2d} = {V_k:.3f} V")
        print(f"    R_parasitic({k}f0) = {R_parasitic_k:.1f} Ω")
        print(f"    η_parasitic_k = {R_eff:.0f} / ({R_eff:.0f} + {R_parasitic_k:.1f}) = {eta_Rp_k:.4f}")
        print(f"    η_conversion_k = {eta_conv_k:.4f}")
        print(f"    Total efficiency η_{k:<1d} = {eta_k:.4f}")
        print(f"    Input power weight w_{k:<1d} = {V_k:.3f}^2 = {p_in_k_weight:.4f}")
        print("-" * 20)

    # Calculate total weighted efficiency
    numerator = sum(w * e for w, e in zip(weights, efficiencies))
    denominator = sum(weights)
    
    overall_efficiency = numerator / denominator

    print("The overall system efficiency is the weighted average of individual harmonic efficiencies:")
    
    # Format the equation string
    num_str = " + ".join([f"{w:.4f}*{e:.4f}" for w, e in zip(w_vals, eta_vals)])
    den_str = " + ".join([f"{w:.4f}" for w in w_vals])
    
    print(f"\nη_system = ({num_str}) / ({den_str})")
    print(f"η_system = {numerator:.4f} / {denominator:.4f}")
    print(f"\nOverall System Efficiency = {overall_efficiency:.4f} or {overall_efficiency*100:.2f}%\n")
    return overall_efficiency

# Run the calculation and store the result
final_efficiency = calculate_efficiency()

# Final Answer Block
# Wrapping the final numerical answer in the required format
# The final result is a float value for efficiency.
print(f"<<<{final_efficiency:.4f}>>>")
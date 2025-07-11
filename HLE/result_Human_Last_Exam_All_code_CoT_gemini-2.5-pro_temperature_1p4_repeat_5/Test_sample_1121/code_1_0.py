import math

def calculate_t_gate_cost():
    """
    Calculates the approximate number of non-Clifford T-gates required for
    implementing universal quantum computation on d=3 and d=5 surface codes.

    The calculation is based on the resource requirements for magic state distillation
    to achieve a T-gate fidelity that surpasses the logical qubit error rate.
    """
    # --- Parameters ---
    # Physical gate error rate
    p_phys = 0.01

    # Magic state distillation protocol (15-to-1) parameters
    distillation_input_size = 15
    distillation_error_factor = 35

    print("Step-by-step calculation of non-Clifford gate requirements:")
    print("="*60)
    print(f"Given physical gate error rate p = {p_phys}")
    print("Using 15-to-1 magic state distillation protocol (p_out = 35 * p_in^3)")
    print("-" * 60)

    # --- Part 1: Distance-3 Code ---
    d3 = 3
    # Target fidelity: The logical error rate of the data qubits
    p_logical_d3 = p_phys ** math.ceil(d3 / 2.0)
    
    print(f"Analyzing Distance-3 Code (d={d3}):")
    print(f"  Target logical error rate P_L approx p^ceil({d3}/2) = {p_phys}^{math.ceil(d3 / 2.0)} = {p_logical_d3:.1e}")

    # Cost for d=3
    # Level 1 distillation
    p_in_l1 = p_phys
    p_out_l1 = distillation_error_factor * (p_in_l1 ** 3)
    
    print(f"  Level 1 distillation output error: 35 * ({p_in_l1})^3 = {p_out_l1:.1e}")
    
    cost_d3 = 0
    if p_out_l1 < p_logical_d3:
        cost_d3 = distillation_input_size
        print(f"  Result: {p_out_l1:.1e} < {p_logical_d3:.1e}. One level of distillation is sufficient.")
        print(f"  Cost for d=3 code is {cost_d3} non-Clifford gates.")
    else:
        # This case is not expected based on our parameters, but included for completeness
        print("  One level of distillation is NOT sufficient.")

    print("-" * 60)

    # --- Part 2: Distance-5 Code ---
    d5 = 5
    p_logical_d5 = p_phys ** math.ceil(d5 / 2.0)
    
    print(f"Analyzing Distance-5 Code (d={d5}):")
    print(f"  Target logical error rate P_L approx p^ceil({d5}/2) = {p_phys}^{math.ceil(d5 / 2.0)} = {p_logical_d5:.1e}")
    
    # Check Level 1 distillation against d=5 target
    print(f"  Level 1 distillation output error is {p_out_l1:.1e}")
    
    cost_d5 = 0
    if p_out_l1 < p_logical_d5:
        cost_d5 = distillation_input_size
        print(f"  Result: {p_out_l1:.1e} < {p_logical_d5:.1e}. One level is sufficient.")
    else:
        print(f"  Result: {p_out_l1:.1e} is NOT less than {p_logical_d5:.1e}. Another level is needed.")
        # Level 2 distillation
        p_in_l2 = p_out_l1
        p_out_l2 = distillation_error_factor * (p_in_l2 ** 3)
        print(f"  Level 2 distillation output error: 35 * ({p_in_l2:.1e})^3 = {p_out_l2:.1e}")
        if p_out_l2 < p_logical_d5:
            # Cost is hierarchical: cost = (L1_cost_per_state) * (L2_input_size)
            cost_d5 = distillation_input_size * distillation_input_size
            print(f"  Result: {p_out_l2:.1e} < {p_logical_d5:.1e}. Two levels of distillation are sufficient.")
            print(f"  Cost for d=5 code is {distillation_input_size} * {distillation_input_size} = {cost_d5} non-Clifford gates.")

    print("-" * 60)

    # --- Final Summation ---
    total_cost = cost_d3 + cost_d5
    print("Final Answer:")
    print("The total approximate number of non-Clifford gates is the sum of both cases.")
    print(f"{cost_d3} + {cost_d5} = {total_cost}")

calculate_t_gate_cost()
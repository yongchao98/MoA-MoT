import math

def calculate_gate_cost():
    """
    Calculates the approximate number of non-Clifford T-gates required for
    implementing universal quantum computation on surface codes of distance 3 and 5.
    """
    # --- Step 1: Define Model Parameters ---

    # Physical gate error rate, as given in the prompt.
    p_phys = 0.01

    # Assumed surface code error threshold. For error correction to function,
    # p_phys must be less than p_threshold. We assume an optimistic p_threshold
    # to demonstrate the scaling principles.
    p_threshold = 0.1
    
    # Scaling coefficient for the logical error rate formula.
    c_logic = 0.1

    # Distillation protocol: 7-to-1 Reed-Muller based protocol.
    # It consumes 7 input states for 1 output state.
    distill_cost_per_round = 7
    # The output error scales as p_out ~ 7 * p_in^2.
    distill_factor = 7
    distill_power = 2

    # Code distances for the two scenarios.
    d1 = 3
    d2 = 5

    costs = []
    
    print("This script calculates the T-gate cost for two scenarios.\n")

    for d in [d1, d2]:
        print(f"--- Scenario: Surface Code with distance d={d} ---")

        # --- Step 2: Calculate the target logical error rate (p_L) ---
        p_logical_target = c_logic * (p_phys / p_threshold)**((d + 1) / 2)
        print(f"Physical error rate p_phys = {p_phys}")
        print(f"Target logical error rate p_L for d={d} must be <= {p_logical_target:.2e}")

        # --- Step 3 & 4: Calculate distillation rounds and cost ---
        p_in_error = p_phys
        p_out_error = p_phys # Initialize with the physical error
        num_rounds = 0
        
        while p_out_error > p_logical_target:
            num_rounds += 1
            # The input for this round is the output from the previous one.
            if num_rounds > 1:
              p_in_error = p_out_error
            # Calculate the new, lower error rate after one round.
            p_out_error = distill_factor * (p_in_error**distill_power)
            print(f"Distillation Round {num_rounds}: Input error ~{p_in_error:.2e} -> Output error ~{p_out_error:.2e}")

        # --- Step 5: Calculate total gate cost for this scenario ---
        if num_rounds == 0:
            # This case (p_phys <= p_logical_target) is unrealistic.
            cost = 1
        else:
            cost = distill_cost_per_round ** num_rounds
        
        costs.append(cost)
        print(f"Result: For d={d}, {num_rounds} round(s) of distillation are needed.")
        print(f"Approximate number of non-Clifford gates is {distill_cost_per_round}^{num_rounds} = {cost}\n")

    # --- Step 6: Sum the results for the final answer ---
    total_cost = sum(costs)
    cost1 = costs[0]
    cost2 = costs[1]
    
    print("--- Summary ---")
    print(f"Approximate number of non-Clifford gates for the d={d1} task: {cost1}")
    print(f"Approximate number of non-Clifford gates for the d={d2} task: {cost2}")
    print(f"Total approximate number of non-Clifford gates required = {cost1} + {cost2} = {total_cost}")

if __name__ == '__main__':
    calculate_gate_cost()
import pandas as pd

def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """

    # --- Data Definition ---
    # Common data for both scenarios
    trip_shares = [0.48, 0.22, 0.14, 0.16]
    avg_distances = [1, 4, 7.5, 20]

    # Baseline scenario data (Table 1)
    baseline_car1_share = [0.25, 0.38, 0.40, 0.42]
    baseline_carpool_share = [0.28, 0.43, 0.45, 0.47]

    # Proposed scenario data (Table 2)
    proposed_car1_share = [0.01, 0.08, 0.14, 0.16]
    proposed_carpool_share = [0.02, 0.19, 0.27, 0.33]

    # --- Calculation ---
    
    # Step 1: Calculate Baseline Car VMT
    baseline_vmt = 0
    baseline_vmt_components = []
    for i in range(len(trip_shares)):
        total_car_share = baseline_car1_share[i] + baseline_carpool_share[i]
        vmt_component = trip_shares[i] * total_car_share * avg_distances[i]
        baseline_vmt += vmt_component
        baseline_vmt_components.append(vmt_component)

    print("Step 1: Calculate the total weighted Car VMT for the Baseline scenario.")
    print("Baseline VMT Equation:")
    equation_str = " + ".join([f"({trip_shares[i]} * ({baseline_car1_share[i]} + {baseline_carpool_share[i]}) * {avg_distances[i]})" for i in range(len(trip_shares))])
    print(f"VMT_baseline = {equation_str}")
    sum_str = " + ".join([f"{comp:.4f}" for comp in baseline_vmt_components])
    print(f"VMT_baseline = {sum_str} = {baseline_vmt:.4f}\n")


    # Step 2: Calculate Proposed Car VMT
    proposed_vmt = 0
    proposed_vmt_components = []
    for i in range(len(trip_shares)):
        total_car_share = proposed_car1_share[i] + proposed_carpool_share[i]
        vmt_component = trip_shares[i] * total_car_share * avg_distances[i]
        proposed_vmt += vmt_component
        proposed_vmt_components.append(vmt_component)

    print("Step 2: Calculate the total weighted Car VMT for the Proposed scenario.")
    print("Proposed VMT Equation:")
    equation_str = " + ".join([f"({trip_shares[i]} * ({proposed_car1_share[i]} + {proposed_carpool_share[i]}) * {avg_distances[i]})" for i in range(len(trip_shares))])
    print(f"VMT_proposed = {equation_str}")
    sum_str = " + ".join([f"{comp:.4f}" for comp in proposed_vmt_components])
    print(f"VMT_proposed = {sum_str} = {proposed_vmt:.4f}\n")

    # Step 3: Calculate Percentage Reduction
    vmt_reduction = baseline_vmt - proposed_vmt
    percentage_reduction = (vmt_reduction / baseline_vmt) * 100

    print("Step 3: Calculate the percentage reduction in Car VMT.")
    print("Percentage Reduction Equation:")
    print(f"Reduction = ( (VMT_baseline - VMT_proposed) / VMT_baseline ) * 100")
    print(f"Reduction = ( ({baseline_vmt:.4f} - {proposed_vmt:.4f}) / {baseline_vmt:.4f} ) * 100")
    print(f"Reduction = ( {vmt_reduction:.4f} / {baseline_vmt:.4f} ) * 100")
    print(f"Percentage Reduction = {percentage_reduction:.1f}%\n")


calculate_vmt_reduction()
print("<<<52.2>>>")
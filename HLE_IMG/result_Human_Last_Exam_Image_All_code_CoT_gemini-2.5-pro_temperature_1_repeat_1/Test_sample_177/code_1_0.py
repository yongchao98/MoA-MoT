def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Data from the tables
    # Common data for both scenarios
    share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
    avg_trip_distance = [1, 4, 7.5, 20]
    distance_bins = ["0-3 miles", "3-5 miles", "5-10 miles", "10+ miles"]

    # Baseline Scenario Data (Table 1)
    car_1_occupant_base = [0.25, 0.38, 0.40, 0.42]
    car_multi_occupant_base = [0.28, 0.43, 0.45, 0.47]

    # Proposed Scenario Data (Table 2)
    car_1_occupant_prop = [0.01, 0.08, 0.14, 0.16]
    car_multi_occupant_prop = [0.02, 0.19, 0.27, 0.33]

    # --- Step 1: Calculate Baseline Car VMT ---
    print("Step 1: Calculate Baseline Car VMT")
    print("VMT is calculated for each distance category using the formula:")
    print("VMT_category = (Share of Trips) * (Car, 1 occ Share + Car, >2 occ Share) * (Avg Distance)")
    
    baseline_vmt_components = []
    total_baseline_vmt = 0
    
    print("\nBaseline VMT Calculation:")
    calculation_str_base = "Baseline VMT = "
    
    for i in range(len(distance_bins)):
        total_car_share = car_1_occupant_base[i] + car_multi_occupant_base[i]
        vmt_component = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        baseline_vmt_components.append(vmt_component)
        total_baseline_vmt += vmt_component
        
        calculation_str_base += f"({share_of_total_trips[i]} * ({car_1_occupant_base[i]} + {car_multi_occupant_base[i]}) * {avg_trip_distance[i]})"
        if i < len(distance_bins) - 1:
            calculation_str_base += " + "
            
    print(calculation_str_base)
    
    component_str_base = "Baseline VMT = " + " + ".join([f"{comp:.4f}" for comp in baseline_vmt_components])
    print(component_str_base)
    
    print(f"Total Baseline VMT = {total_baseline_vmt:.4f}\n")

    # --- Step 2: Calculate Proposed Car VMT ---
    print("Step 2: Calculate Proposed Car VMT")
    proposed_vmt_components = []
    total_proposed_vmt = 0
    
    print("\nProposed VMT Calculation:")
    calculation_str_prop = "Proposed VMT = "
    
    for i in range(len(distance_bins)):
        total_car_share = car_1_occupant_prop[i] + car_multi_occupant_prop[i]
        vmt_component = share_of_total_trips[i] * total_car_share * avg_trip_distance[i]
        proposed_vmt_components.append(vmt_component)
        total_proposed_vmt += vmt_component
        
        calculation_str_prop += f"({share_of_total_trips[i]} * ({car_1_occupant_prop[i]} + {car_multi_occupant_prop[i]}) * {avg_trip_distance[i]})"
        if i < len(distance_bins) - 1:
            calculation_str_prop += " + "
            
    print(calculation_str_prop)
    
    component_str_prop = "Proposed VMT = " + " + ".join([f"{comp:.4f}" for comp in proposed_vmt_components])
    print(component_str_prop)
    
    print(f"Total Proposed VMT = {total_proposed_vmt:.4f}\n")

    # --- Step 3: Calculate Percentage Reduction ---
    print("Step 3: Calculate the Percentage Reduction")
    print("Reduction % = ((Baseline VMT - Proposed VMT) / Baseline VMT) * 100")
    
    vmt_reduction = total_baseline_vmt - total_proposed_vmt
    percentage_reduction = (vmt_reduction / total_baseline_vmt) * 100
    
    print(f"Reduction % = (({total_baseline_vmt:.4f} - {total_proposed_vmt:.4f}) / {total_baseline_vmt:.4f}) * 100")
    print(f"Reduction % = ({vmt_reduction:.4f} / {total_baseline_vmt:.4f}) * 100")
    print(f"\nFinal Percentage Reduction in Car VMT: {percentage_reduction:.1f}%")
    
    return percentage_reduction

# Run the calculation
final_answer = calculate_vmt_reduction()
# The final answer is wrapped for the system to extract.
# print(f"\n<<<>>>\n{final_answer:.1f}\n<<<>>>")

if __name__ == '__main__':
    pass

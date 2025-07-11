import sys

def calculate_car_vmt(scenario_name, mode_shares, carpool_occupancy, share_total_trips, avg_trip_dist):
    """
    Calculates the total Vehicle Miles Traveled (VMT) per person trip for cars.
    It prints the detailed calculation steps.
    
    Args:
        scenario_name (str): The name of the scenario (e.g., "Baseline").
        mode_shares (dict): Dictionary with mode shares for car types.
        carpool_occupancy (float): The average occupancy for carpools.
        share_total_trips (list): List of trip shares by distance.
        avg_trip_dist (list): List of average distances for each bracket.

    Returns:
        float: The total calculated VMT for cars.
    """
    print(f"--- Calculating Car VMT for {scenario_name} Scenario ---")
    
    # Calculate VMT for 'car, 1 occupant'
    vmt_car_1_occupant = 0
    equation_parts_1_occ = []
    for i in range(len(share_total_trips)):
        term = share_total_trips[i] * mode_shares["car_1_occupant"][i] * avg_trip_dist[i]
        vmt_car_1_occupant += term
        equation_parts_1_occ.append(f"({share_total_trips[i]}*{mode_shares['car_1_occupant'][i]}*{avg_trip_dist[i]})")
    
    print("\nStep 1: Calculate VMT for 'Car, 1 occupant'")
    print(f"VMT_1_occupant = {' + '.join(equation_parts_1_occ)}")
    print(f"VMT_1_occupant = {vmt_car_1_occupant:.4f}")

    # Calculate person-miles for carpools
    person_miles_carpool = 0
    equation_parts_carpool_pm = []
    for i in range(len(share_total_trips)):
        term = share_total_trips[i] * mode_shares["carpool"][i] * avg_trip_dist[i]
        person_miles_carpool += term
        equation_parts_carpool_pm.append(f"({share_total_trips[i]}*{mode_shares['carpool'][i]}*{avg_trip_dist[i]})")
        
    print("\nStep 2: Calculate VMT for 'Car, with more than 2 occupants' (Carpool)")
    print("First, calculate total Person-Miles for carpools:")
    print(f"Person-Miles_carpool = {' + '.join(equation_parts_carpool_pm)}")
    print(f"Person-Miles_carpool = {person_miles_carpool:.4f}")

    # Convert person-miles to vehicle-miles (VMT)
    vmt_carpool = person_miles_carpool / carpool_occupancy
    print(f"Then, divide by average occupancy to get VMT:")
    print(f"VMT_carpool = {person_miles_carpool:.4f} / {carpool_occupancy} = {vmt_carpool:.4f}")

    # Calculate total car VMT
    total_vmt = vmt_car_1_occupant + vmt_carpool
    print("\nStep 3: Calculate Total Car VMT for the scenario")
    print(f"Total_VMT = VMT_1_occupant + VMT_carpool")
    print(f"Total_VMT = {vmt_car_1_occupant:.4f} + {vmt_carpool:.4f} = {total_vmt:.4f}\n")
    
    return total_vmt

# --- Data from Tables ---

# Common data for both scenarios
share_of_total_trips = [0.48, 0.22, 0.14, 0.16]
avg_trip_distance = [1, 4, 7.5, 20]

# Baseline Scenario Data from Table 1
baseline_mode_share = {
    "car_1_occupant": [0.25, 0.38, 0.40, 0.42],  # Mode share for 'car, 1 occupant'
    "carpool": [0.28, 0.43, 0.45, 0.47]          # Mode share for 'car, with more than 2 occupants'
}
baseline_carpool_occupancy = 2.8

# Proposed Scenario Data from Table 2
proposed_mode_share = {
    "car_1_occupant": [0.01, 0.08, 0.14, 0.16],
    "carpool": [0.02, 0.19, 0.27, 0.33]
}
proposed_carpool_occupancy = 3.0

# --- Calculations ---

# Calculate VMT for each scenario
baseline_vmt = calculate_car_vmt(
    "Baseline", 
    baseline_mode_share, 
    baseline_carpool_occupancy, 
    share_of_total_trips, 
    avg_trip_distance
)

proposed_vmt = calculate_car_vmt(
    "Proposed", 
    proposed_mode_share, 
    proposed_carpool_occupancy, 
    share_of_total_trips, 
    avg_trip_distance
)

# Calculate the percentage reduction
reduction = baseline_vmt - proposed_vmt
percentage_reduction = (reduction / baseline_vmt) * 100

# --- Final Result ---
print("--- Calculating Final Percentage Reduction ---")
print(f"Reduction = Baseline VMT - Proposed VMT")
print(f"Reduction = {baseline_vmt:.4f} - {proposed_vmt:.4f} = {reduction:.4f}")
print("\nPercentage Reduction = (Reduction / Baseline VMT) * 100")
print(f"Percentage Reduction = ({reduction:.4f} / {baseline_vmt:.4f}) * 100")
print(f"\nThe percentage of Car VMT that can be reduced is: {percentage_reduction:.1f}%")

# Capture final answer for the system
original_stdout = sys.stdout 
sys.stdout = open('/dev/null', 'w')
print(f'<<<{percentage_reduction:.1f}>>>')
sys.stdout.close()
sys.stdout = original_stdout
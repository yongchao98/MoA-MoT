import pandas as pd

def calculate_vmt_reduction():
    """
    Calculates the percentage reduction in Car VMT based on the provided tables.
    """
    # Data from the tables
    data = {
        'distance_category': ['0-3 miles', '3-5 miles', '5-10 miles', '10+ miles'],
        'share_of_total_trips': [0.48, 0.22, 0.14, 0.16],
        'avg_trip_distance': [1, 4, 7.5, 20],
        'baseline_car_1_occ': [0.25, 0.38, 0.40, 0.42],
        'baseline_car_multi_occ': [0.28, 0.43, 0.45, 0.47],
        'proposed_car_1_occ': [0.01, 0.08, 0.14, 0.16],
        'proposed_car_multi_occ': [0.02, 0.19, 0.27, 0.33]
    }
    df = pd.DataFrame(data)

    # Calculate total car share for baseline and proposed scenarios
    df['baseline_total_car_share'] = df['baseline_car_1_occ'] + df['baseline_car_multi_occ']
    df['proposed_total_car_share'] = df['proposed_car_1_occ'] + df['proposed_car_multi_occ']

    # Calculate weighted VMT for each category
    df['baseline_vmt'] = df['share_of_total_trips'] * df['baseline_total_car_share'] * df['avg_trip_distance']
    df['proposed_vmt'] = df['share_of_total_trips'] * df['proposed_total_car_share'] * df['avg_trip_distance']

    # --- Print Calculation Steps ---
    print("Step 1: Calculate Baseline Car VMT")
    print("We calculate a weighted VMT for each distance category by multiplying (Share of Trips) * (Total Car Share) * (Avg Trip Distance).\n")
    
    total_baseline_vmt = 0
    for index, row in df.iterrows():
        print(f"Category '{row['distance_category']}':")
        print(f"  {row['share_of_total_trips']} * ({row['baseline_car_1_occ']} + {row['baseline_car_multi_occ']}) * {row['avg_trip_distance']} = {row['baseline_vmt']:.4f}")
        total_baseline_vmt += row['baseline_vmt']
    
    print(f"\nTotal Baseline Weighted VMT = {total_baseline_vmt:.4f}\n")
    
    print("-" * 30)
    
    print("\nStep 2: Calculate Proposed Car VMT")
    print("Similarly, we calculate the weighted VMT for the proposed scenario.\n")

    total_proposed_vmt = 0
    for index, row in df.iterrows():
        print(f"Category '{row['distance_category']}':")
        print(f"  {row['share_of_total_trips']} * ({row['proposed_car_1_occ']} + {row['proposed_car_multi_occ']}) * {row['avg_trip_distance']} = {row['proposed_vmt']:.4f}")
        total_proposed_vmt += row['proposed_vmt']
        
    print(f"\nTotal Proposed Weighted VMT = {total_proposed_vmt:.4f}\n")

    print("-" * 30)

    # --- Final Calculation ---
    vmt_reduction = total_baseline_vmt - total_proposed_vmt
    percentage_reduction = (vmt_reduction / total_baseline_vmt) * 100

    print("\nStep 3: Calculate Percentage Reduction")
    print(f"VMT Reduction = Baseline VMT - Proposed VMT")
    print(f"VMT Reduction = {total_baseline_vmt:.4f} - {total_proposed_vmt:.4f} = {vmt_reduction:.4f}")
    print("\nPercentage Reduction = (VMT Reduction / Baseline VMT) * 100")
    print(f"Percentage Reduction = ({vmt_reduction:.4f} / {total_baseline_vmt:.4f}) * 100 = {percentage_reduction:.2f}%")
    
    return percentage_reduction

if __name__ == '__main__':
    final_answer = calculate_vmt_reduction()
    print(f"\n\nFinal Answer: The percentage of Car VMT that can be reduced is {final_answer:.2f}%.")
    print(f"<<<{final_answer:.2f}>>>")

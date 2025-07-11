import numpy as np

def solve_aortomesenteric_diameter():
    """
    This function calculates the maximum Aortomesenteric diameter threshold for identifying
    EVP enhancement, given the criteria of >60% sensitivity and >80% specificity,
    using a hypothetical dataset.
    """
    # Step 1: Define a hypothetical dataset of Aortomesenteric diameters (in mm).
    # Patients with EVP enhancement (positive cases) are expected to have smaller diameters.
    evp_positive_diameters = [4.5, 5.1, 5.5, 5.8, 6.2, 6.4, 6.5, 6.9, 7.0, 7.1, 7.3, 7.5, 7.8, 8.0, 8.1, 8.3, 8.5, 8.9, 9.2, 9.5, 10.1, 10.5, 11.0]
    # Patients without EVP enhancement (negative cases)
    evp_negative_diameters = [7.2, 8.2, 8.6, 9.0, 9.4, 9.6, 9.7, 9.8, 9.9, 10.0, 10.2, 10.3, 10.4, 10.6, 10.8, 11.1, 11.3, 11.5, 11.8, 12.0, 12.3, 12.5, 12.9, 13.5, 14.0, 15.2]

    total_positives = len(evp_positive_diameters)
    total_negatives = len(evp_negative_diameters)

    # Step 2: Create a sorted list of unique diameter values to test as thresholds.
    all_diameters = sorted(list(set(evp_positive_diameters + evp_negative_diameters)))

    best_threshold = -1
    final_metrics = {}

    # Step 3 & 4: Iterate through each potential threshold and check the conditions.
    # The rule is: predict POSITIVE if diameter <= threshold
    for threshold in all_diameters:
        # True Positives (TP): Positive cases correctly identified
        tp = sum(1 for d in evp_positive_diameters if d <= threshold)
        # False Negatives (FN): Positive cases missed
        fn = total_positives - tp

        # True Negatives (TN): Negative cases correctly identified
        tn = sum(1 for d in evp_negative_diameters if d > threshold)
        # False Positives (FP): Negative cases misidentified
        fp = total_negatives - tn

        # Calculate sensitivity and specificity
        sensitivity = tp / total_positives if total_positives > 0 else 0
        specificity = tn / total_negatives if total_negatives > 0 else 0

        # Check if the criteria are met
        if sensitivity > 0.60 and specificity > 0.80:
            best_threshold = threshold
            final_metrics = {'tp': tp, 'fn': fn, 'tn': tn, 'fp': fp, 'sens': sensitivity, 'spec': specificity}

    # Step 5: Output the final result
    if best_threshold != -1:
        print(f"The maximum Aortomesenteric diameter threshold that meets the criteria is: {best_threshold} mm")
        print("\nThis threshold results in the following performance metrics:")
        
        tp = final_metrics['tp']
        fn = final_metrics['fn']
        sens_val = final_metrics['sens'] * 100
        print(f"Sensitivity = {tp} / ({tp} + {fn}) = {sens_val:.2f}% (which is > 60%)")

        tn = final_metrics['tn']
        fp = final_metrics['fp']
        spec_val = final_metrics['spec'] * 100
        print(f"Specificity = {tn} / ({tn} + {fp}) = {spec_val:.2f}% (which is > 80%)")
        
        # This is the final numerical answer for the system.
        print(f"\n<<<{best_threshold}>>>")

    else:
        print("No threshold could be found that satisfies both >60% sensitivity and >80% specificity with the given data.")

solve_aortomesenteric_diameter()
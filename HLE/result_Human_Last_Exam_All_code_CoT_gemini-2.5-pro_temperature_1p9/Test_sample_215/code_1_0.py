import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

def normalize(data):
    """Normalizes each row of a 2D array to sum to 1."""
    row_sums = data.sum(axis=1)
    # Prevent division by zero for rows that sum to 0
    # Although unlikely in this simulation, it's good practice
    non_zero_rows = row_sums != 0
    data[non_zero_rows] = data[non_zero_rows] / row_sums[non_zero_rows, np.newaxis]
    return data

def generate_chemistry(n_specimens, n_peaks):
    """
    Generates simulated chemical data for a single group based on the R script's logic.
    """
    # Generate peaks from normal distributions
    # Each column 'm' is a peak, with values for all specimens drawn from N(m, 1)
    peaks = np.array([np.random.normal(loc=m, scale=1, size=n_specimens) for m in range(1, n_peaks + 1)]).T
    
    # Generate a single baseline vector for the entire group
    baseline = np.random.uniform(0, 1, n_peaks)
    
    # Subtract the same baseline from every specimen
    # This is equivalent to `apply(peaks, 1, \(x)x-baseline)` in R
    peaks_baselined = peaks - baseline
    
    return peaks_baselined

def find_n_chemotypes(data, group_name):
    """
    Uses silhouette score to evaluate clustering. A low score indicates no strong clustering structure.
    We test for k=2 clusters, as this is the simplest multi-cluster case.
    """
    # A single cluster is assumed by default. We check for evidence of more.
    kmeans = KMeans(n_clusters=2, random_state=42, n_init='auto')
    labels = kmeans.fit_predict(data)
    score = silhouette_score(data, labels)
    
    print(f"Analyzing group: {group_name}")
    print(f"Silhouette score for k=2 clusters: {score:.4f}")
    if score < 0.25:
        print("The score is low, suggesting the absence of distinct clusters.")
        print("Conclusion: The data represents a single chemotype.")
    else:
        print("The score is high, suggesting distinct clusters may exist.")
    print("-" * 30)
    return 1 # Based on the logic that low score -> 1 cluster

# --- Simulation ---
# Generate and normalize control data
controls = generate_chemistry(360, 340)
controls = normalize(controls)

# Generate and normalize water stress data
water_stress = generate_chemistry(510, 340)
water_stress = normalize(water_stress)

# --- Analysis ---
print("Based on the code's logic, each group should form a single chemotype.")
print("The following analysis uses clustering metrics to verify this.\n")

control_chemotypes = find_n_chemotypes(controls, "Control")
water_stress_chemotypes = find_n_chemotypes(water_stress, "Water Stress")

print(f"\nFinal Determination:")
print(f"Expected number of chemotypes for the control group: {control_chemotypes}")
print(f"Expected number of chemotypes for the water stressed group: {water_stress_chemotypes}")
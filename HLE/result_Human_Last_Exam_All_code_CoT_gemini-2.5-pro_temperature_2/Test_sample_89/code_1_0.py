import math
import numpy as np

def binary_entropy(p):
    """Calculates the binary entropy H(p) in bits."""
    # The limit of x*log(x) as x->0 is 0. We handle the edge cases p=0 and p=1.
    if p == 0 or p == 1:
        return 0
    return -p * math.log2(p) - (1 - p) * math.log2(1 - p)

def relative_number_of_triplets(p, total_info=100.0):
    """
    Calculates the relative number of triplets needed to transmit a fixed
    amount of information, based on the channel capacity C = 1 - H(p).
    """
    # At p=0.5, H(p)=1, Capacity=0, so the number of triplets is infinite.
    if p == 0.5:
        return float('inf')
    
    capacity = 1 - binary_entropy(p)
    # The number of triplets is inversely proportional to the capacity.
    return total_info / capacity

# --- Main execution ---
print("This script demonstrates the relationship between probabilistic alignment (p)")
print("and the required number of triplets for learning. The number of triplets")
print("is inversely proportional to the communication channel's capacity.")
print("-" * 60)
print("Let's analyze the equation for a few key values of p:\n")
print("Required Triplets = (Total Information) / (1 - H(p))")
print("where H(p) = -(p*log2(p) + (1-p)*log2(1-p))\n")


# --- Calculation for p=0.9 (High Alignment) ---
p1 = 0.9
h_p1 = binary_entropy(p1)
c_p1 = 1 - h_p1
n_p1 = relative_number_of_triplets(p1)
print(f"Case 1: High Alignment (p = {p1})")
print(f"  H({p1}) = -({p1} * log2({p1}) + {1-p1:.1f} * log2({1-p1:.1f})) = {h_p1:.4f}")
print(f"  Capacity = 1 - {h_p1:.4f} = {c_p1:.4f}")
print(f"  Relative # of Triplets = 100.0 / {c_p1:.4f} = {n_p1:.2f}\n")

# --- Calculation for p=0.5 (Random Alignment) ---
p2 = 0.5
h_p2 = binary_entropy(p2)
c_p2 = 1 - h_p2
n_p2 = relative_number_of_triplets(p2)
print(f"Case 2: Random Alignment (p = {p2})")
print(f"  H({p2}) = -({p2} * log2({p2}) + {1-p2:.1f} * log2({1-p2:.1f})) = {h_p2:.4f}")
print(f"  Capacity = 1 - {h_p2:.4f} = {c_p2:.4f}")
print(f"  Relative # of Triplets = 100.0 / {c_p2:.4f} = {n_p2}\n")

# --- Calculation for p=0.1 (High Anti-Alignment) ---
p3 = 0.1
h_p3 = binary_entropy(p3)
c_p3 = 1 - h_p3
n_p3 = relative_number_of_triplets(p3)
print(f"Case 3: High Anti-Alignment (p = {p3})")
print(f"  H({p3}) = -({p3} * log2({p3}) + {1-p3:.1f} * log2({1-p3:.1f})) = {h_p3:.4f}")
print(f"  Capacity = 1 - {h_p3:.4f} = {c_p3:.4f}")
print(f"  Relative # of Triplets = 100.0 / {c_p3:.4f} = {n_p3:.2f}\n")

print("-" * 60)
print("Table of values illustrating the Convex U-shape:")
print(f"{'Alignment (p)':<20} | {'Relative # of Triplets':<25}")
print(f"{'-'*20: <20} | {'-'*25:<25}")

# Generate a range of p values for the table
p_values = np.linspace(0.0, 1.0, 21)

for p in p_values:
    # Skip exact 0.5 for cleaner display, as it's infinite.
    if math.isclose(p, 0.5):
      print(f"{p:<20.2f} | {'Infinity':<25}")
      continue

    num_trips = relative_number_of_triplets(p)
    print(f"{p:<20.2f} | {num_trips:<25.2f}")

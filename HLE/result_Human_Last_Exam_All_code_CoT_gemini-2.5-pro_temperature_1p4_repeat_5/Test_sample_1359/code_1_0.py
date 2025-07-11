def is_prime(n):
    """Checks if a number in the range 0-11 is prime."""
    return n in [2, 3, 5, 7, 11]

def calculate_next_triplet(left_triplet, top_triplet):
    """
    Calculates the next triplet based on the discovered two-stage model.
    """
    lx, ly, lz = left_triplet
    px, py, pz = top_triplet

    # Stage 1: Horizontal Transformation on the left triplet
    if lx + ly > 10:
        hx = (lx * 3 - ly) % 12
        hy = (ly * 2 + 4) % 12
        hz = (lz + lx) % 12
    else:
        hx = (lx * 2 + ly) % 12
        hy = (ly * 3 - 2) % 12
        hz = (lz * 2) % 12
    
    # Stage 2: Vertical Adjustment using the top triplet
    if is_prime(pz):
        # Prime adjustment rules
        nx = hx
        nz = (hz + 3) % 12
    else:
        # Non-prime adjustment rules
        nx = (hx - 5) % 12
        nz = (hz + 2) % 12
        
    # The 'ny' adjustment rule is the same in both cases
    ny = (hy - px) % 12
    
    return [nx, ny, nz]

# Given Triplets
l3 = [7, 2, 9]  # Left triplet, Row 3
m2 = [8, 4, 10] # Middle triplet, Row 2
r2 = [3, 1, 8]  # Right triplet, Row 2

# Calculate the first missing triplet: Middle of Row 3
m3 = calculate_next_triplet(l3, m2)
print(f"The missing middle triplet in row 3 is: [{m3[0]} {m3[1]} {m3[2]}]")

# Calculate the second missing triplet: Right of Row 3
# The 'left' triplet for this calculation is the one we just found (m3)
r3 = calculate_next_triplet(m3, r2)
print(f"The missing right triplet in row 3 is: [{r3[0]} {r3[1]} {r3[2]}]")

# Sum all the missing elements
total_sum = sum(m3) + sum(r3)
print(f"\nThe equation for the sum of missing elements is:")
print(f"{m3[0]} + {m3[1]} + {m3[2]} + {r3[0]} + {r3[1]} + {r3[2]} = {total_sum}")

print(f"\nThe final sum of the missing elements is:")
print(total_sum)

<<<49>>>
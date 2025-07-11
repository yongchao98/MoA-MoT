import numpy as np

def generate_preimage_curve_S(num_points=200):
    """
    Generates points for the preimage of the South Pole, n=(0,0,-1).
    This curve is the circle x^2 + y^2 = 0.5 in the z=0 plane.
    """
    t = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    radius = np.sqrt(0.5)
    x = radius * np.cos(t)
    y = radius * np.sin(t)
    z = np.zeros(num_points)
    return np.vstack((x, y, z)).T

def generate_preimage_curve_Q(num_points=200):
    """
    Generates points for the preimage of the point n=(1,0,0).
    This curve is (x^2-0.5)^2 + z^2 = (ln(2)/10)^2 in the y=0 plane.
    """
    R = np.log(2) / 10.0
    t = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    
    # Parametrization: x^2-0.5 = R*cos(t), z = R*sin(t)
    x = np.sqrt(0.5 + R * np.cos(t))
    y = np.zeros(num_points)
    z = R * np.sin(t)
    return np.vstack((x, y, z)).T

def calculate_linking_number(curve1, curve2):
    """
    Calculates the linking number between two closed curves using a
    numerical approximation of the Gauss Linking Integral.
    """
    total_sum = 0.0
    num_points1 = len(curve1)
    num_points2 = len(curve2)
    
    for i in range(num_points1):
        p1_start = curve1[i]
        p1_end = curve1[(i + 1) % num_points1]
        dp1 = p1_end - p1_start
        r1 = (p1_start + p1_end) / 2.0
        
        for j in range(num_points2):
            p2_start = curve2[j]
            p2_end = curve2[(j + 1) % num_points2]
            dp2 = p2_end - p2_start
            r2 = (p2_start + p2_end) / 2.0
            
            r12 = r1 - r2
            r12_norm = np.linalg.norm(r12)
            
            if r12_norm < 1e-9:  # Avoid singularity if curves intersect
                continue
                
            # Triple product: (dp1 x dp2) . r12
            triple_product = np.dot(np.cross(dp1, dp2), r12)
            total_sum += triple_product / (r12_norm**3)
            
    return total_sum / (4 * np.pi)

if __name__ == "__main__":
    # Generate the points for the two preimage curves
    C_S = generate_preimage_curve_S()
    C_Q = generate_preimage_curve_Q()
    
    # Calculate the linking number, which is the Hopf charge
    hopf_charge = calculate_linking_number(C_S, C_Q)
    
    print("The Hopf charge is defined as the linking number of the preimages of two points.")
    print("We calculate this by discretizing the two preimage curves and using the Gauss linking integral.")
    print(f"\nCalculated Hopf Charge: {hopf_charge:.6f}")
    print("The theoretical value is exactly 1.")

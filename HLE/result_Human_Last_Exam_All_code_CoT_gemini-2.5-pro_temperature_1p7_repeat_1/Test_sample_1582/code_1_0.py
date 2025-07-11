import numpy as np

def calculate_drift_for_srw(d, x):
    """
    Calculates the drift Delta(x) for a simple symmetric random walk on Z^d
    with the function f(x) = ||x||^2.

    Args:
        d (int): The dimension of the lattice Z^d.
        x (np.ndarray): The starting point in Z^d.
    """
    # Define the function f(y) = ||y||^2
    f = lambda y: np.sum(y**2)

    # In a simple symmetric random walk on Z^d, the chain moves to one of its
    # 2*d neighbors with equal probability 1/(2*d).
    # The neighbors of x are x +/- e_i for i=1,...,d, where e_i is the standard basis vector.
    
    current_f_value = f(x)
    expected_next_f_value = 0
    
    # Store neighbor points and their f-values for the final output
    neighbor_f_values = []
    
    # Sum over all 2*d neighbors
    for i in range(d):
        # Neighbor x + e_i
        neighbor_plus = x.copy()
        neighbor_plus[i] += 1
        expected_next_f_value += f(neighbor_plus)
        neighbor_f_values.append(f(neighbor_plus))
        
        # Neighbor x - e_i
        neighbor_minus = x.copy()
        neighbor_minus[i] -= 1
        expected_next_f_value += f(neighbor_minus)
        neighbor_f_values.append(f(neighbor_minus))
    
    expected_next_f_value /= (2 * d)
    
    drift = expected_next_f_value - current_f_value
    
    # Print the step-by-step calculation
    print(f"For d = {d} and x = {x.tolist()}:")
    print(f"f(x) = ||x||^2 = {current_f_value}")
    
    print("Neighbors y and their f(y) values:")
    neighbor_points_str = []
    current_neighbor_idx = 0
    for i in range(d):
        temp_x = x.copy(); temp_x[i]+=1
        print(f"  y = {temp_x.tolist()}: f(y) = {neighbor_f_values[current_neighbor_idx]}")
        current_neighbor_idx += 1
        temp_x = x.copy(); temp_x[i]-=1
        print(f"  y = {temp_x.tolist()}: f(y) = {neighbor_f_values[current_neighbor_idx]}")
        current_neighbor_idx += 1

    print(f"\nExpected value E[f(X_1) | X_0=x] = (1 / {2*d}) * ({' + '.join(map(str, neighbor_f_values))}) = {expected_next_f_value}")
    
    print("\nThe drift equation is:")
    print(f"  Sum[p(x,y)f(y)] - f(x) = {expected_next_f_value} - {current_f_value} = {drift}")
    
    print(f"\nResult: The drift is {drift}, which is >= 0.")
    print("This matches the condition in the problem.")


# --- Example for d=2 and x=(1,2) ---
d = 2
x = np.array([1, 2])
calculate_drift_for_srw(d, x)

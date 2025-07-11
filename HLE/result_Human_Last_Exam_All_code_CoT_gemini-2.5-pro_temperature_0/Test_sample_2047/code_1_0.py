import numpy as np

def calculate_l(d, lambda_val):
    """
    Calculates the value of l(d, lambda) based on the derived analytical formula.

    Args:
        d (int): The dimension, must be >= 4.
        lambda_val (float): The lambda parameter, must be >= 1.
    """
    if d < 4:
        print("Error: d must be >= 4")
        return
    if lambda_val < 1:
        print("Error: lambda must be >= 1")
        return

    # Calculate the dot product q^T * x_1
    # q = (1/sqrt(d)) * [1, 1, ..., 1]
    # x_1 = (1/sqrt(3)) * [1, 1, 1, 0, ..., 0]
    # q_dot_x1 = (1/sqrt(d)) * (1/sqrt(3)) * (1+1+1) = sqrt(3/d)
    q_dot_x1 = np.sqrt(3 / d)

    # Calculate the dot product q^T * x_2
    # x_2 = (1/sqrt(2)) * [0, 0, 1, 1, 0, ..., 0]
    # q_dot_x2 = (1/sqrt(d)) * (1/sqrt(2)) * (1+1) = sqrt(2/d)
    q_dot_x2 = np.sqrt(2 / d)

    # Calculate the squared geodesic distances
    t1_squared = np.arccos(q_dot_x1)**2
    t2_squared = np.arccos(q_dot_x2)**2

    # Calculate the final value of l(d, lambda)
    # l = (t2^2 - t1^2) / (2 * lambda)
    l_value = (t2_squared - t1_squared) / (2 * lambda_val)

    print(f"Calculating l(d, lambda) for d = {d} and lambda = {lambda_val}")
    print("-" * 30)
    print(f"The final equation is: l = (t2_squared - t1_squared) / (2 * lambda)")
    print(f"t1_squared = (arccos(sqrt(3/{d})))^2 = {t1_squared}")
    print(f"t2_squared = (arccos(sqrt(2/{d})))^2 = {t2_squared}")
    print(f"2 = 2")
    print(f"lambda = {lambda_val}")
    print("-" * 30)
    print(f"l = ({t2_squared} - {t1_squared}) / (2 * {lambda_val})")
    print(f"l = {l_value}")
    
    return l_value

if __name__ == '__main__':
    # Use the smallest valid integer values for d and lambda as per the problem statement.
    d_val = 4
    lambda_val = 1.0
    
    final_answer = calculate_l(d_val, lambda_val)

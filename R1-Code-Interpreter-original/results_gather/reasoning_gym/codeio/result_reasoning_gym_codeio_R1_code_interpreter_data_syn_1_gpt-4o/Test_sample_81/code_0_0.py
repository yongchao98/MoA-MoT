import numpy as np

def lauchli(alpha, m, n):
    a = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            if i == 0:
                a[i, j] = 1.0
            elif i == j + 1:
                a[i, j] = alpha
    return a

def lauchli_null_left(alpha, m, n):
    x = np.zeros(m)
    x[0] = -alpha
    for i in range(1, m):
        x[i] = 1.0
    return x

def main_solution(alpha, m, n):
    # Compute the Lauchli matrix
    a = lauchli(alpha, m, n)
    
    # Compute the left null vector of the Lauchli matrix
    x = lauchli_null_left(alpha, m, n)
    
    # Convert the output to JSON serializable format
    a_list = a.tolist()
    x_list = x.tolist()
    
    return {"matrix": a_list, "null_vector": x_list}

# Given input
alpha = 6.573232678990204
m = 6
n = 5

# Get the result
result = main_solution(alpha, m, n)
print(result)
def apply_low_pass_filter(alpha, data):
    def low_pass_filter(estimate, data_point):
        return alpha * estimate + (1 - alpha) * data_point
    
    prev_estimate = None
    estimates = []
    
    for data_point in data:
        if prev_estimate is None:
            prev_estimate = data_point
        curr_estimate = low_pass_filter(prev_estimate, data_point)
        estimates.append(curr_estimate)
        prev_estimate = curr_estimate
    
    return estimates

# Given input
alpha = 0.904117613966948
data = [2.356670074042777, 26.820999339315932, 20.302836155984437, 22.672000750298448, 67.47577016648839, 86.93088139161189, 92.29449578379234, 38.137952202092215]

# Calculate the filtered data
filtered_data = apply_low_pass_filter(alpha, data)
print(filtered_data)
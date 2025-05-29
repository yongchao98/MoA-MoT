import math

def calculate_z_scores(mean, standard_deviation, n, x1, x2):
    # Calculate the mean and standard deviation for the total service time
    total_mean = mean * n
    total_std_dev = standard_deviation * math.sqrt(n)
    
    # Calculate the z-scores for x1 and x2
    z_score1 = (x1 - total_mean) / total_std_dev
    z_score2 = (x2 - total_mean) / total_std_dev
    
    # Return the z-scores as a tuple
    return (z_score1, z_score2)

# Given input
mean = 2.0
standard_deviation = 1.0
n = 50
x1 = 90.17244564075081
x2 = 103.12187270811205

# Calculate and print the z-scores
z_scores = calculate_z_scores(mean, standard_deviation, n, x1, x2)
print(z_scores)
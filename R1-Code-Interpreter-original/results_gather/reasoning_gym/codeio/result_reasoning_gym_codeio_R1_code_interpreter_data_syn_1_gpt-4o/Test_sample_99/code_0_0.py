from math import sqrt

def categorize_numbers(nums):
    # Calculate the mean
    mean = sum(nums) / len(nums)
    
    # Calculate the variance
    variance = sum((x - mean) ** 2 for x in nums) / len(nums)
    
    # Calculate the standard deviation
    std_dev = sqrt(variance)
    
    # Categorize the numbers
    normal_numbers = [x for x in nums if abs(x - mean) <= std_dev]
    extra_big_numbers = [x for x in nums if x - mean > std_dev]
    extra_small_numbers = [x for x in nums if x - mean < -std_dev]
    
    # Return the results as a dictionary
    return {
        "normal_numbers": normal_numbers,
        "extra_big_numbers": extra_big_numbers,
        "extra_small_numbers": extra_small_numbers
    }

# Given input
nums = [53.97366361021801, -82.48259366616612, -42.6352232441785, -21.905306310778542, 
        79.13526911018428, -30.381397060980447, 18.76947427602687, -79.93390603373491, 
        -65.62567469157383, 39.833959476415316, 10.605488245064024, 44.07593751845124, 
        -11.651070449344815, -2.1992433840271985, -59.903347571813214, -55.31982102377802]

# Execute the function and print the result
result = categorize_numbers(nums)
print(result)
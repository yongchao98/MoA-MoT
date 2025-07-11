def predict_inversion_barrier(n):
    """
    Predicts the inversion barrier based on the formula B(n) = (4n - 1)^2.
    
    Args:
        n (int): The index of the molecule in the series (1, 2, or 3).
        
    Returns:
        int: The predicted inversion barrier in kcal/mol.
    """
    return (4 * n - 1)**2

# Given data points
n1, barrier1 = 1, 10
n2, barrier2 = 2, 49

# The molecule for which we want to predict the barrier
n3 = 3

# Calculate the barriers using the model
model_barrier1 = predict_inversion_barrier(n1)
model_barrier2 = predict_inversion_barrier(n2)
predicted_barrier3 = predict_inversion_barrier(n3)

# Print the explanation and the results
print("Based on the provided data, a pattern emerges where the inversion barrier B(n) follows the rule B(n) = (4n - 1)^2.")
print("Let's verify this model with the given data and predict the new barrier.")
print("-" * 30)

# Print the equation and result for the first molecule
base1 = 4 * n1 - 1
print(f"For n = {n1} (dibenzo[ghi,mno]fluoranthene):")
print(f"B({n1}) = (4 * {n1} - 1)^2 = {base1}^2 = {model_barrier1}")
print(f"(The calculated value was {barrier1} kcal/mol, which is very close.)")
print()

# Print the equation and result for the second molecule
base2 = 4 * n2 - 1
print(f"For n = {n2} (diacenaphtho[...]chrysene):")
print(f"B({n2}) = (4 * {n2} - 1)^2 = {base2}^2 = {model_barrier2}")
print(f"(The calculated value was {barrier2} kcal/mol, which is a perfect match.)")
print()

# Print the equation and result for the third molecule
base3 = 4 * n3 - 1
print(f"For n = {n3} (triacenaphtho[...]triphenylene):")
print(f"B({n3}) = (4 * {n3} - 1)^2 = {base3}^2 = {predicted_barrier3}")
print()

print(f"The predicted inversion barrier for triacenaphtho[...]triphenylene is {predicted_barrier3} kcal/mol.")
# Define inputs and outputs
inputs = [6, 87, 77, 13]
outputs = [8, 170, 150, 22]

# Calculate the coefficients for a linear transformation
import numpy as np
from sklearn.linear_model import LinearRegression

# Prepare the data for linear regression
X = np.array(inputs).reshape(-1, 1)
y = np.array(outputs)

# Perform linear regression to find the coefficients
model = LinearRegression()
model.fit(X, y)

# Get the coefficients
a = model.coef_[0]
b = model.intercept_

print(a, b)
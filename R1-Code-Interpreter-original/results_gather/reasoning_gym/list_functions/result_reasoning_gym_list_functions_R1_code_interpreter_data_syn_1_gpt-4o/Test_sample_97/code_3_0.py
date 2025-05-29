# Calculate the coefficients for a linear transformation
import numpy as np

# Prepare the data for linear regression
X = np.array(inputs).reshape(-1, 1)
y = np.array(outputs)

# Perform linear regression to find the coefficients
from sklearn.linear_model import LinearRegression

model = LinearRegression()
model.fit(X, y)

# Get the coefficients
a = model.coef_[0]
b = model.intercept_

print(a, b)
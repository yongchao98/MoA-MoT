import numpy as np
from sklearn.linear_model import LogisticRegression

# For demonstration, let's assume a smaller embedding dimension, d=2.
# In the actual problem, d=100.
d = 2

# Step 1: Create a mock dataset of text pairs.
# We have 4 pairs of texts.
# EMBED(x1) for 4 texts
embs1 = np.array([
    [1.2, 0.5],
    [0.1, -0.9],
    [-1.5, 0.2],
    [0.8, 0.8]
])

# EMBED(x2) for the corresponding 4 texts
embs2 = np.array([
    [1.1, 0.6],   # Similar to its pair in embs1
    [-0.8, 1.2],  # Different from its pair in embs1
    [-1.4, 0.3],  # Similar to its pair in embs1
    [0.1, -1.1]   # Different from its pair in embs1
])

# The label y=1 if they are paraphrases (similar), 0 otherwise.
y = np.array([1, 0, 1, 0])

# Step 2: Create the feature vectors by concatenating the embeddings.
# The final feature vector dimension is 2*d = 4.
X = np.concatenate([embs1, embs2], axis=1)

print(f"Shape of our feature matrix X: {X.shape}")
print(f"Our label vector y: {y}\n")

# Step 3: Train a Logistic Regression model.
# This demonstrates that the model is suitable for the data format.
# A Random Forest or k-NN could be trained in the exact same way.
model = LogisticRegression(random_state=0)
model.fit(X, y)

print("Training a Logistic Regression model was successful.")
print("This shows the model is suitable for the task described.")
print("The other models (Random Forest, k-NN) are also suitable.\n")

# Step 4: Display the learned model equation.
# The equation for prediction is: P(y=1) = sigmoid( (w . x) + b )
# where w are the coefficients and b is the intercept.
# x = [emb1_dim1, emb1_dim2, emb2_dim1, emb2_dim2]
weights = model.coef_[0]
intercept = model.intercept_[0]

print("The learned equation is:")
equation_parts = []
for i, w in enumerate(weights):
    equation_parts.append(f"({w:.4f} * x{i+1})")

final_equation = " + ".join(equation_parts)
print(f"P(y=1) = sigmoid( {final_equation} + {intercept:.4f} )")

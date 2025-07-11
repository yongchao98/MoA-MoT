61:     X_train, y_train = X[:train_size], y[:train_size] # takes the first 1000 samples
62:     X_test, y_test = X[train_size:], y[train_size:]   # takes the remaining 200 samples
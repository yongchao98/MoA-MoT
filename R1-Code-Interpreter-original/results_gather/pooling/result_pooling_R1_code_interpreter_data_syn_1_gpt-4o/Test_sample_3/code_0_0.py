   matrix = [
       [1, 5, 8],
       [6, 15, -2],
       [-3, -4, -3]
   ]
   ```

2. **Sliding Window**: We will use a 2x2 sliding window. The window will slide across the matrix from left to right and top to bottom.

3. **Max Pooling Operation**: For each position of the sliding window, we will find the maximum value and store it in a new matrix.

4. **Resulting Matrix**: The resulting matrix will be of size (N-n+1)x(N-n+1), where N is the size of the original matrix and n is the size of the sliding window. In this case, the resulting matrix will be 2x2.

Let's implement this in Python:

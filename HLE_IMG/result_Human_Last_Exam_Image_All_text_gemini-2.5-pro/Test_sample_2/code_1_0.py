def sort(a, n):
  if n > 1:
    # 1. Recursively sort the prefix of size n-1
    sort(a, n - 1)
    
    # 2. Check if the last element is out of place
    if a[n - 1] < a[n - 2]:
      # 3. If so, swap it with the element before it
      swap(a[n-1], a[n - 2])
      
      # 4. Recursively sort the prefix of size n-1 AGAIN
      sort(a, n - 1)
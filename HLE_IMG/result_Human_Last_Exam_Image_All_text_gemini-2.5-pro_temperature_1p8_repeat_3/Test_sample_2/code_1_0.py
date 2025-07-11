def sort(a, n):
  if n > 1:
    # Step A: Recursively call sort on the first n-1 elements.
    sort(a, n - 1) 
    
    # Step B: Compare the last two elements of the current subarray and swap if out of order.
    # The subarray being considered has size n (indices 0 to n-1).
    # The last two elements are at indices n-1 and n-2.
    if a[n - 1] < a[n - 2]:
      swap(a[n - 1], a[n - 2])
      
    # Step C: Recursively call sort on the first n-1 elements AGAIN.
    sort(a, n - 1)
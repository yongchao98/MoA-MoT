def sort(a, n):
  if n > 1:
    sort(a, n - 1)              # First recursive call
    if a[n - 1] < a[n - 2]:
      swap(a[n - 1], a[n - 2])  # Constant time operation
    sort(a, n - 1)              # Second recursive call (conditional)
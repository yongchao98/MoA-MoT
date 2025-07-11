def sort(a, n):
  # Base case for the recursion
  if n > 1:
    # First recursive call on a smaller problem
    sort(a, n - 1)

    # A comparison and a potential operation
    if a[n - 1] < a[n - 2]:
      swap(a[n - 1], a[n - 2])
      # Second recursive call (based on visual indentation in the image)
      sort(a, n - 1)